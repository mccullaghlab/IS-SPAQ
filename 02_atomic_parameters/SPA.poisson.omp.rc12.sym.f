CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM fit_gr
        include 'omp_lib.h'
        integer nmax,mpts,hmax,mtyp
        parameter(nmax=200,mpts=8000000,hmax=120,mtyp=100)
        double precision y(mpts),x(3,mpts),wy(mpts)
        double precision rmin(mtyp)
        double precision xpos(3,nmax)
        double precision g(mtyp*hmax)
        double precision g2(mtyp*hmax),g3(mtyp*hmax)
        integer ni(mtyp*hmax)
        integer imap(mtyp*hmax)
        double precision nc(mtyp*hmax,mtyp*hmax)
        double precision nc2(mtyp*hmax,mtyp*hmax)
        integer ityp(nmax),inow(nmax)
        double precision xnow,y0,e0,s0,p1,p2,p3,dx,rcut2
        integer i,j,npts,ntyp,npar,k,jp,kp,nmat,igo,ncpu,jt
        character*64 fname
C
        call omp_set_dynamic(.false.)
        ncpu=6
        dx=0.1d0
        rcut2=12.d0**2
C read in the molecule coordinates and types
        write(fname,888) 
888     format('ADI.crd')
        open(20,FILE=fname,STATUS='old')
        read(20,*) npar,ntyp
        do i=1,npar
          read(20,*) ityp(i),xpos(1,i),xpos(2,i),xpos(3,i)
        enddo
        close(20)
C
        y0=0.018962962963d0
        e0=dlog(y0)-0.577215665d0
        s0=3.1415926535898d0**2/6.d0
C
        do i=1,ntyp*hmax
          do j=1,ntyp*hmax
            nc(j,i)=0.d0
          enddo
          g(i)=0.d0
          ni(i)=0
        enddo
C read in the distribution function
        write(fname,887)
887     format('ADI.p2.d0.09.gr3')
        open(20,FILE=fname,STATUS='old')
        read(20,*) npts
        do i=1,npts
          read(20,*) x(1,i),x(2,i),x(3,i),y(i),p1,p2,p3
          jp=int(y(i)/y0+0.1d0)
          y(i)=e0
          wy(i)=s0
          do j=1,jp
            y(i)=y(i)+1.d0/dble(j)
            wy(i)=wy(i)-1.d0/dble(j*j)
          enddo
          wy(i)=1.d0/wy(i)
          do j=1,npar
            jt=ityp(j)
            xnow=(x(1,i)-xpos(1,j))**2+
     &           (x(2,i)-xpos(2,j))**2+
     &           (x(3,i)-xpos(3,j))**2
            if(xnow.lt.rcut2) then
              k=int(dsqrt(xnow)/dx)+1
              k=hmax*(jt-1)+k
              ni(k)=ni(k)+int(wy(i)+0.1d0)
            endif
          enddo
        enddo
        close(20)
C
        write(6,*) 'finding rmin'
        do i=1,ntyp
          jp=(i-1)*hmax+1
          j=1
          do while(ni(jp).le.1.and.j.lt.hmax)
            jp=jp+1
            j=j+1
          enddo
          if(j.eq.hmax) then
            rmin(i)=1.d10
          else
            rmin(i)=dble(j-1)*dx
            rmin(i)=rmin(i)*rmin(i)
          endif
        enddo
C
        write(0,*) 'looping over grid'
        do i=1,npts
          if(wy(i).gt.0.1d0) then
            igo=1
            j=1
            jt=ityp(j)
            do while (igo.gt.0)
              xnow=(x(1,i)-xpos(1,j))**2
     &             +(x(2,i)-xpos(2,j))**2
     &             +(x(3,i)-xpos(3,j))**2
              if(xnow.gt.rmin(jt)) then
                xnow=dsqrt(xnow)
                inow(j)=int(xnow/dx)+1
                if(inow(j).le.hmax) then
                  inow(j)=hmax*(jt-1)+inow(j)
                  igo=igo+1
                else
                  inow(j)=0
                endif
              else
                igo=0
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
              do j=1,npar
                jp=inow(j)
                if(jp.gt.0) then
                  g(jp)=g(jp)+wy(i)*y(i)
                  do k=1,npar
                    kp=inow(k)
                    if(kp.gt.0) then
                      nc(kp,jp)=nc(kp,jp)+wy(i)
                    endif
                  enddo
                endif
              enddo
            endif
          endif
        enddo
C
        write(0,*) 'reordering matrix'
        nmat=0
        do i=1,hmax*ntyp
          if(nc(i,i).gt.0.5d0) then
            nmat=nmat+1
            if(nmat.ne.i) then
              do j=1,hmax*ntyp
                nc(j,nmat)=nc(j,i)
              enddo
              do j=1,hmax*ntyp
                nc(nmat,j)=nc(i,j)
              enddo
              g(nmat)=g(i)
            endif
            imap(i)=nmat
          else
            imap(i)=0
          endif
        enddo
C     
        write(0,*) 'total:',nmat
        do i=1,nmat
          do j=1,nmat
            nc2(j,i)=nc(j,i)
          enddo
          g2(i)=g(i)
        enddo
        call inv(nc,g,nmat,ncpu)
C
        do i=1,nmat
          g3(i)=-g2(i)
          do j=1,nmat
            g3(i)=g3(i)+nc2(j,i)*g(j)
          enddo
        enddo
C
        write(fname,886) 
886     format('sample.SPA.dat')
        open(20,FILE=fname,STATUS='unknown')
        jp=0
        do i=1,ntyp
          do j=1,hmax
            jp=jp+1
            xnow=(j-0.5d0)*dx
            kp=imap(jp)
            if(kp.eq.0) then
              write(20,*) xnow,0.d0,i,0.d0,-1
            else
              write(20,*) xnow,g(kp),i,g3(kp),ni(jp)
            endif
          enddo
        enddo
        close(20)
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
