CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM fit_gr
        include 'omp_lib.h'
        integer nmax,hmax,nhist,mtyp,mcpu
        parameter(nmax=200,hmax=120,nhist=3000000,mtyp=100,mcpu=20)
        double precision y,wy
        double precision q(nmax)
        double precision xg(3,nhist),ng(nhist)
        double precision x(3,nmax)
        double precision dr(3,nmax),edotr(nmax)
        double precision e(mtyp*hmax)
        double precision a(mtyp*hmax,mtyp*hmax)
        double precision b(mtyp*hmax),b0(mtyp*hmax)
        double precision ao(mtyp*hmax,mtyp*hmax,mcpu)
        double precision bo(mtyp*hmax,mcpu)
        integer ityp(nmax),jtyp(nmax),inow(nmax),imin(2,mtyp)
        double precision xnow,xnow2,p1,p2,p3,x1,x2,x3
        double precision dnow,dmax,gnow,r0,r2,rdotr
        double precision dx1,dx2,dx3,ei,pi,p2i,enow(3)
        double precision rho,eps,p0,dx,rcut2
        integer ix,ip,ir,istop,mpts,tid,n,m
        integer i,j,npts,ntyp,npar,k,jp,kp,nmat,igo,ncpu,jt
        character*64 fname
C
        call omp_set_dynamic(.false.)
        ncpu=6
        dx=0.1d0
        rcut2=12.d0**2
C
        rho=0.0074d0
        p0=0.2303d0
        eps=2.3473d0
C read in the molecule coordinates and types
        write(fname,888)
888     format('ADI.crd')
        open(20,FILE=fname,STATUS='old')
        read(20,*) npar,ntyp
        do i=1,npar
          read(20,*) ityp(i),x(1,i),x(2,i),x(3,i)
          jtyp(ityp(i))=i
        enddo
        close(20)
C
        open(20,FILE='q.dat',STATUS='old')
        do i=1,npar
          read(20,*) q(i)
C          q(i)=q(i)*3.d0/4.d0/3.1415926536d0/0.0075d0/0.2915d0
C     &         /dsqrt(332.d0)*(1.d0-1.d0/4.39d0)
          q(i)=q(i)*3.d0/4.d0/3.1415926536d0/rho/p0
     &         /dsqrt(332.d0)*(1.d0-1.d0/eps)
        enddo
        close(20)
C read in previous fit - determine imin
        write(fname,884)
884     format('sample.SPA2.dat')
        open(20,FILE=fname,STATUS='unknown')
        nmat=0
        do i=1,ntyp
          imin(1,i)=nmat
          do j=1,hmax
            read(20,*) xnow,gnow,jt,m
            if(m.lt.0) then
              imin(2,i)=j
            elseif(m.ne.0) then
              nmat=nmat+1
              e(nmat)=-0.001d0*q(jtyp(i))/xnow**2
            else
              write(0,*) 'There is a gap in the fit.'
            endif
          enddo
        enddo
        close(20)
C
        do i=1,nmat
          b0(i)=0.d0
        enddo
C read in the distribution function
        mpts=0
        write(6,*) 'reading p(r)'
        write(0,*) 'reading p(r)'
        write(fname,887)
887     format('sample.gr3')
        open(20,FILE=fname,STATUS='old')
        read(20,*) npts
        do i=1,npts
          read(20,*) x1,x2,x3,y,p1,p2,p3,(pi,j=1,5),wy
C          read(20,*) x1,x2,x3,y,p1,p2,p3,wy
          if(wy.gt.0.5d0) then
            igo=0
            do j=1,npar
              jt=ityp(j)
              dx1=x1-x(1,j)
              dx2=x2-x(2,j)
              dx3=x3-x(3,j)
              r2=dx1**2+dx2**2+dx3**2
              if (r2.lt.rcut2) then
                r0=dsqrt(r2)
                ir=int(r0/dx)+1
                if(ir.gt.imin(2,jt)) then
                  k=imin(1,jt)+ir-imin(2,jt)
                  edotr(j)=(p1*dx1+p2*dx2+p3*dx3)/r0*wy
                  inow(j)=k
                  igo=igo+1
                else
                  igo=-npar
                endif
              else
                inow(j)=0
              endif
            enddo
            if(igo.gt.0) then
              mpts=mpts+1
              xg(1,mpts)=x1
              xg(2,mpts)=x2
              xg(3,mpts)=x3
              ng(mpts)=wy
              do j=1,npar
                k=inow(j)
                if(k.gt.0) then
                  b0(k)=b0(k)-edotr(j)
                endif
              enddo
            endif
          endif
        enddo
        close(20)
C
        istop=0
        do while (istop.eq.0)
C
          write(6,*) 'looping over grid: ', mpts
          write(0,*) 'looping over grid: ', mpts
          
C$OMP parallel num_threads(ncpu) default(none)
C$OMP& private(dr,edotr,inow,enow,tid,ix,i,j,k,jt,x1,x2,x3,wy,
C$OMP&             xnow2,xnow,ip,jp,kp,ei,pi,p2i,rdotr,n)
C$OMP& shared(xg,ng,ao,bo,a,b,b0,x,q,ityp,imin,e,mpts,npar,nmat,ncpu,dx)
          tid=omp_get_thread_num()+1
          do i=1,nmat
            do j=1,nmat
              ao(j,i,tid)=0.d0
            enddo
            bo(i,tid)=0.d0
          enddo
C$OMP do schedule(guided)
          do ix=1,mpts
            x1=xg(1,ix)
            x2=xg(2,ix)
            x3=xg(3,ix)
            wy=ng(ix)
C
            enow(1)=0.d0
            enow(2)=0.d0
            enow(3)=0.d0
            do j=1,npar
              jt=ityp(j)
              dr(1,j)=x1-x(1,j)
              dr(2,j)=x2-x(2,j)
              dr(3,j)=x3-x(3,j)
              xnow2=dr(1,j)**2+dr(2,j)**2+dr(3,j)**2
              xnow=dsqrt(xnow2)
              dr(1,j)=dr(1,j)/xnow
              dr(2,j)=dr(2,j)/xnow
              dr(3,j)=dr(3,j)/xnow
              ip=int(xnow/dx)+1
              if(ip.le.hmax) then
                ip=imin(1,jt)+ip-imin(2,jt)
                enow(1)=enow(1)+dr(1,j)*e(ip)
                enow(2)=enow(2)+dr(2,j)*e(ip)
                enow(3)=enow(3)+dr(3,j)*e(ip)
                inow(j)=ip
              else
C                enow(1)=enow(1)-dr(1,j)*q(j)/xnow2
C                enow(2)=enow(2)-dr(2,j)*q(j)/xnow2
C                enow(3)=enow(3)-dr(3,j)*q(j)/xnow2
                inow(j)=0
              endif
              enow(1)=enow(1)-dr(1,j)*q(j)/xnow2
              enow(2)=enow(2)-dr(2,j)*q(j)/xnow2
              enow(3)=enow(3)-dr(3,j)*q(j)/xnow2
            enddo
C
            ei=enow(1)**2+enow(2)**2+enow(3)**2
            ei=dsqrt(ei)
            enow(1)=enow(1)/ei
            enow(2)=enow(2)/ei
            enow(3)=enow(3)/ei
            pi=1.d0/dtanh(ei)-1.d0/ei
            p2i=1.d0-pi**2-3.d0*pi/ei
C
            do j=1,npar
              edotr(j)=enow(1)*dr(1,j)+enow(2)*dr(2,j)+enow(3)*dr(3,j)
            enddo
C
            do j=1,npar
              jp=inow(j)
              if(jp.gt.0) then
                bo(jp,tid)=bo(jp,tid)+edotr(j)*pi*wy
                do k=1,npar
                  kp=inow(k)
                  if(kp.gt.0) then
                    rdotr=dr(1,j)*dr(1,k)
     &                   +dr(2,j)*dr(2,k)
     &                   +dr(3,j)*dr(3,k)
                    ao(kp,jp,tid)=ao(kp,jp,tid)
     &                   +(edotr(j)*edotr(k)*p2i+rdotr*pi/ei)*wy
                  endif
                enddo
              endif
            enddo
          enddo
C$OMP enddo
C$OMP do schedule(guided)
          do i=1,nmat
            b(i)=b0(i)
            do n=1,ncpu
              b(i)=b(i)+bo(i,n)
            enddo
            do j=1,nmat
              a(j,i)=0.d0
              do n=1,ncpu
                a(j,i)=a(j,i)+ao(j,i,n)
              enddo
            enddo
          enddo
C$OMP enddo
C$OMP end parallel
C
          write(6,*) 'total:',nmat
          write(0,*) 'total:',nmat
          call inv(a,b,nmat,ncpu)
          dmax=0.d0
          do i=1,nmat
            dnow=dabs(b(i))
            if(dnow.gt.dmax) dmax=dnow
C            e(i)=e(i)-b(i)
          enddo
          if(dmax.lt.5.d0) then
            do i=1,nmat
              e(i)=e(i)-b(i)
            enddo
          else
            do i=1,nmat
              e(i)=e(i)-b(i)/dmax*5.d0
            enddo
          endif
C
          write(6,*) dmax
          write(0,*) dmax
          if (dmax.lt.0.01d0) istop=1
C
          write(fname,886)
886       format('sample.pol2.dat')
          open(20,FILE=fname,STATUS='unknown')
          jp=0
          do i=1,ntyp
            jt=jtyp(i)
            do j=1,imin(2,i)
              xnow=(j-0.5d0)*dx
              write(20,*) xnow,0.d0,i,0.d0
            enddo
            do j=imin(2,i)+1,hmax
              jp=jp+1
              xnow=(j-0.5d0)*dx
              write(20,*) xnow,e(jp),i,q(jt)
            enddo
          enddo
          close(20)
        enddo
C     
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine inv(mat,val,nmat,ncpu)
        include 'omp_lib.h'
        integer nmax,hmax
        parameter(nmax=21,hmax=120)
        double precision mat(nmax*hmax,nmax*hmax)
        double precision val(nmax*hmax)
        integer irow2(20),icol2(20)
        double precision big2(20)
        integer ipiv(nmax*hmax)
        integer i,j,k,irow,icol,nmat,ncpu,tid
        integer indxc(nmax*hmax),indxr(nmax*hmax)
        double precision big,pivinv,dum
C
        do i=1,nmat
          ipiv(i)=0
        enddo
C
        do i=1,nmat
C          write(0,*) i
          do j=1,ncpu
            big2(j)=0.d0
          enddo
C$OMP parallel num_threads(ncpu) default(none)
C$OMP& private(j,k,tid)
C$OMP& shared(i,nmat,ipiv,mat,big2,irow2,icol2)
          tid=omp_get_thread_num()+1
C$OMP do schedule(guided)
          do j=1,nmat
            if(ipiv(j).eq.0) then
              do k=1,nmat
                if(ipiv(k).eq.0) then
                  if(abs(mat(j,k)).ge.big2(tid)) then
                    big2(tid)=abs(mat(j,k))
                    irow2(tid)=j
                    icol2(tid)=k
                  endif
                endif
              enddo
            endif
          enddo
C$OMP enddo          
C$OMP end parallel
          big=big2(1)
          irow=irow2(1)
          icol=icol2(1)
          do j=2,ncpu
            if(big2(j).gt.big) then
              big=big2(j)
              irow=irow2(j)
              icol=icol2(j)
            endif
          enddo
          ipiv(icol)=ipiv(icol)+1
C
          if(irow.ne.icol) then
C$OMP parallel do schedule(static) num_threads(ncpu) default(none)
C$OMP& private(j,dum)
C$OMP& shared(nmat,mat,irow,icol)
            do j=1,nmat
              dum=mat(irow,j)
              mat(irow,j)=mat(icol,j)
              mat(icol,j)=dum
            enddo
C$OMP end parallel do
            dum=val(irow)
            val(irow)=val(icol)
            val(icol)=dum
          endif
          indxr(i)=irow
          indxc(i)=icol
          pivinv=1.d0/mat(icol,icol)
          mat(icol,icol)=1.d0
C$OMP parallel do schedule(static) num_threads(ncpu) default(none)
C$OMP& private(j,dum)
C$OMP& shared(nmat,mat,icol,pivinv)
          do j=1,nmat
            mat(icol,j)=mat(icol,j)*pivinv
          enddo
C$OMP end parallel do
          val(icol)=val(icol)*pivinv
C
C$OMP parallel do schedule(guided) num_threads(ncpu) default(none)
C$OMP& private(j,dum,k)
C$OMP& shared(nmat,mat,icol,val)
          do j=1,nmat
            if(j.ne.icol) then
              dum=mat(j,icol)
              mat(j,icol)=0.d0
              do k=1,nmat
                mat(j,k)=mat(j,k)-mat(icol,k)*dum
              enddo
              val(j)=val(j)-val(icol)*dum
            endif
          enddo
C$OMP end parallel do
        enddo
C
C        do i=nmat,1,-1
C          if(indxr(i).ne.indxc(i)) then
C            do j=1,nmat
C              dum=mat(j,indxr(i))
C              mat(j,indxr(i))=mat(j,indxc(i))
C              mat(j,indxc(i))=dum
C            enddo
C          endif
C        enddo
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
