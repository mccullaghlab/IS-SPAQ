CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM fit_gr
        include 'omp_lib.h'
        integer nmax,hmax,mtyp,mcpu
        parameter(nmax=200,hmax=120,mtyp=100,mcpu=20)
        double precision y,hy0,dy
        double precision rmin(mtyp)
        double precision x(3,nmax)
        double precision g(mtyp*hmax)
        double precision a(mtyp*hmax,mtyp*hmax)
        double precision ao(mtyp*hmax,mtyp*hmax,mcpu)
        double precision b(mtyp*hmax),b0(mtyp*hmax)
        double precision bo(mtyp*hmax,mcpu)
        integer ityp(nmax),inow(nmax),imin(2,nmax)
        double precision xnow,p1,p2,p3,x1,x2,x3
        double precision dnow,dmax,gnow,r0,r2,dx,rcut2
        double precision sum
        integer igrd,ix,iy,iz,ip,ir,istop,m,tid,isum
        integer i,j,npts,nx1,nx2,ntyp,npar,k,jp,kp,nmat,igo,ncpu,jt
        character*64 fname
        integer ib(mtyp*hmax),ibo(mtyp*hmax,mcpu)
C
C     the density was messed up a touch...
        dy=1.d0
C        
        call omp_set_dynamic(.false.)
        ncpu=6
        dx=0.1d0
        rcut2=12.d0**2
C  Half the density reported in a cell with one count.
        hy0=0.018939629340d0*0.5d0
C read in the molecule coordinates and types
        write(fname,888) 
888     format('ADI.crd')
        open(20,FILE=fname,STATUS='old')
        read(20,*) npar,ntyp
        do i=1,npar
          read(20,*) ityp(i),x(1,i),x(2,i),x(3,i)
        enddo
        close(20)
C read in previous fit - determine imin
        write(fname,884) 
884     format('sample.SPA.dat')
        open(20,FILE=fname,STATUS='unknown')
        nmat=0
        do i=1,ntyp
          imin(1,i)=nmat
          do j=1,hmax
            read(20,*) xnow,gnow,jt,dnow,m
            if(m.lt.0) then
              imin(2,i)=j
              rmin(i)=xnow
            elseif(m.ne.0) then
              nmat=nmat+1
              g(nmat)=gnow
            else
              write(0,*) 'there is a gap in the fit'
            endif
          enddo
          rmin(i)=rmin(i)+0.05d0
          rmin(i)=rmin(i)*rmin(i)
        enddo
        close(20)
C
        do i=1,nmat
          b0(i)=0.d0
        enddo
C read in the distribution function
        write(0,*) 'reading g(r)'
        write(fname,887)
887     format('ADI.p2.d0.09.gr3')
        open(20,FILE=fname,STATUS='old')
        read(20,*) npts
        do i=1,npts
          read(20,*) x1,x2,x3,y,p1,p2,p3
          if(y.gt.hy0) then
            y=y*dy
            do j=1,npar
              jt=ityp(j)
              r2=(x1-x(1,j))**2+(x2-x(2,j))**2+(x3-x(3,j))**2
              if (r2.lt.rcut2.and.r2.gt.rmin(jt)) then
                r0=dsqrt(r2)
                ir=int(r0/dx)+1
                k=imin(1,jt)+ir-imin(2,jt)
                b0(k)=b0(k)-y
              endif
            enddo
          endif
        enddo
        close(20)
C
        istop=0
        do while (istop.eq.0)
C
          write(0,*) 'looping over grid'
          nx1=int(npts**(1.d0/3.d0)+0.01d0)
          if(nx1**3.ne.npts) write(0,*) 'ERROR: is the gr3 a cube?'
          nx2=nx1**2
          npts=npts-1
C$OMP parallel num_threads(ncpu) default(none)
C$OMP& private(inow,tid,igrd,i,ix,iy,iz,x1,x2,x3,igo,gnow,j,jt,xnow,
C$OMP&          ip,jp,kp,m,sum,isum)
C$OMP& shared(ao,bo,ibo,a,b,b0,ib,g,x,rmin,ityp,imin,npar,nmat,ncpu,
C$OMP&          npts,nx1,nx2,dx,rcut2)
          tid=omp_get_thread_num()+1
          do i=1,nmat
            do j=1,nmat
              ao(j,i,tid)=0.d0
            enddo
            bo(i,tid)=0.d0
            ibo(i,tid)=0
          enddo
C$OMP do schedule(guided)
          do igrd=0,npts
            ix=igrd/nx2+1
            iy=(igrd-(ix-1)*nx2)/nx1+1
            iz=mod(igrd,nx1)+1
            x1=(dble(ix)-0.5d0)*0.25d0-25.d0
            x2=(dble(iy)-0.5d0)*0.25d0-25.d0
            x3=(dble(iz)-0.5d0)*0.25d0-25.d0
C
            igo=1
            gnow=0.d0
            j=1
            jt=ityp(j)
            do while (igo.gt.0)
              xnow=(x(1,j)-x1)**2+(x(2,j)-x2)**2+(x(3,j)-x3)**2
              if(xnow.lt.rcut2) then
                if(xnow.gt.rmin(jt)) then
                  xnow=dsqrt(xnow)
                  ip=int(xnow/dx)+1
                  if(ip.le.hmax) then
                    ip=imin(1,jt)+ip-imin(2,jt)
                    gnow=gnow+g(ip)
                    inow(j)=ip
                    igo=igo+1
                  else
                    inow(j)=0
                  endif
                else
                  igo=0
                endif
              else
                inow(j)=0
              endif
              if(j.eq.npar) then
                if (igo.gt.1) then
                  igo=-1
                else
                  igo=0
                endif
              else
                j=j+1
                jt=ityp(j)
              endif
            enddo
C
            if(igo.eq.-1) then
              gnow=dexp(gnow)
              do j=1,npar
                jp=inow(j)
                if(jp.gt.0) then
                  bo(jp,tid)=bo(jp,tid)+gnow
                  ibo(jp,tid)=ibo(jp,tid)+1
                  do k=1,npar
                    kp=inow(k)
                    if(kp.gt.0) then
                      ao(kp,jp,tid)=ao(kp,jp,tid)+gnow
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
C$OMP enddo
C$OMP do schedule(static)
          do i=1,nmat
            sum=b0(i)
            isum=0
            do m=1,ncpu
              sum=sum+bo(i,m)
              isum=isum+ibo(i,m)
            enddo
            b(i)=sum
            ib(i)=isum
            do j=1,nmat
              sum=0.d0
              do m=1,ncpu
                sum=sum+ao(j,i,m)
              enddo
              a(j,i)=sum
            enddo
          enddo
C$OMP enddo
C$OMP end parallel
C
          write(0,*) 'total:',nmat
          call inv(a,b,nmat,ncpu)
          dmax=0.d0
          do i=1,nmat
            dnow=dabs(b(i))
            if(dnow.gt.dmax) dmax=dnow
            g(i)=g(i)-b(i)
          enddo
C     
          write(0,*) dmax
          if (dmax.lt.0.01d0) istop=1
C
          write(fname,886)
886       format('sample.SPA2.dat')
          open(20,FILE=fname,STATUS='unknown')
          jp=0
          do i=1,ntyp
            do j=1,imin(2,i)
              xnow=(j-0.5d0)*dx
              write(20,*) xnow,0.d0,i,-1
            enddo
            do j=imin(2,i)+1,hmax
              jp=jp+1
              xnow=(j-0.5d0)*dx
              write(20,*) xnow,g(jp),i,ib(jp)
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
