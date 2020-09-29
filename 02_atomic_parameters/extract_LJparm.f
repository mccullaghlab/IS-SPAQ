CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM LJparams
        integer nmax,hmax
        parameter(nmax=30000,hmax=120)
C paramtop values
        integer natom,nbonh,nbona,ntheth,ntheta
        integer nphih,nphia,nnb,nres,ntyp,nmol
        double precision charge(nmax),mass(nmax)
        integer ityp(nmax),nbparm(nmax),nrpnt(nmax),nmpnt(nmax)
        double precision kbnd(nmax),rbnd(nmax),kang(nmax),tang(nmax)
        double precision kdih(nmax),pdih(nmax),ndih(nmax)
        double precision scee(nmax),scnb(nmax)
        double precision alj(nmax),blj(nmax)
        integer bhlist(3,nmax),balist(3,nmax)
        integer ahlist(4,nmax),aalist(4,nmax)
        integer dhlist(5,nmax),dalist(5,nmax)
        integer excl(nmax),nexc(nmax)
        double precision kumb1,kumb2,rumb,mtot1,mtot2
        integer umblist(nmax,2)
C
        character*64 prmtop
        integer i,it,nlj,nsolu
C
        common /param/ charge,mass,kbnd,rbnd,kang,tang,kdih,pdih,ndih,
     &       scee,scnb,alj,blj,kumb1,kumb2,rumb,mtot1,mtot2,ityp,
     &       nbparm,nrpnt,nmpnt,bhlist,balist,ahlist,aalist,dhlist,
     &       dalist,excl,nexc,umblist,natom,ntyp,nbonh,nbona,ntheth,
     &       ntheta,nphih,nphia,nnb,nmol,nres
C
        write(prmtop,999)
999     format("ADI.prmtop")
        call setparam(prmtop)
        nsolu=22
C
        open(20,FILE="LJparm.dat",STATUS="unknown")
        write(6,*) "#    rmin    epsilon - for symmetric pairs"
        write(6,*) "# Solute"
        do i=1,nsolu
          it=ityp(i)
          nlj=ntyp*(it-1)+it
          nlj=nbparm(nlj)
          write(6,888) (alj(nlj)/blj(nlj))**(1.d0/6.d0),
     &         blj(nlj)**2/alj(nlj)/12.d0
        enddo
888     format(f10.6,1X,f10.6)
        write(6,*) "# Solvent"
        do i=nsolu+1,nsolu+5
          it=ityp(i)
          nlj=ntyp*(it-1)+it
          nlj=nbparm(nlj)
          write(6,888) (alj(nlj)/blj(nlj))**(1.d0/6.d0),
     &         blj(nlj)**2/alj(nlj)/12.d0
        enddo
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine setparam(prmtop)
        integer nmax
        parameter(nmax=30000)
        integer natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia,natyp
        integer nnb,nres,numbnd,numang,nptra,nmol
        double precision charge(nmax),mass(nmax)
        integer ityp(nmax),nbparm(nmax),nrpnt(nmax),nmpnt(nmax)
        double precision kbnd(nmax),rbnd(nmax),kang(nmax),tang(nmax)
        double precision kdih(nmax),pdih(nmax),ndih(nmax)
        double precision scee(nmax),scnb(nmax)
        double precision alj(nmax),blj(nmax)
        integer bhlist(3,nmax),balist(3,nmax)
        integer ahlist(4,nmax),aalist(4,nmax)
        integer dhlist(5,nmax),dalist(5,nmax)
        integer excl(nmax),nexc(nmax)
        double precision kumb1,kumb2,rumb,mtot1,mtot2
        integer umblist(nmax,2)
        integer i,j,inext,i1,i2,i3,i4
        character*64 prmtop
C
        common /param/ charge,mass,kbnd,rbnd,kang,tang,kdih,pdih,ndih,
     &       scee,scnb,alj,blj,kumb1,kumb2,rumb,mtot1,mtot2,ityp,
     &       nbparm,nrpnt,nmpnt,bhlist,balist,ahlist,aalist,dhlist,
     &       dalist,excl,nexc,umblist,natom,ntyp,nbonh,nbona,ntheth,
     &       ntheta,nphih,nphia,nnb,nmol,nres
C
        open(20,FILE=prmtop,STATUS='old')
C  numbers
        do i=1,6
          read(20,*)
        enddo
        read(20,999) natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia
        read(20,998) nnb,nres,numbnd,numang,nptra,natyp
C  charge
        inext=7+(natom-1)/20
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/5+1
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) charge(j),charge(j+1),charge(j+2),
     &         charge(j+3),charge(j+4)
        enddo
C  mass
        inext=5+(natom-1)/10
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/5+1
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) mass(j),mass(j+1),mass(j+2),
     &         mass(j+3),mass(j+4)
        enddo
C  ityp
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/10+1
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) ityp(j),ityp(j+1),ityp(j+2),ityp(j+3),ityp(j+4),
     &         ityp(j+5),ityp(j+6),ityp(j+7),ityp(j+8),ityp(j+9)
        enddo
C  n-excluded
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=(natom-1)/10+1
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) nexc(j),nexc(j+1),nexc(j+2),nexc(j+3),nexc(j+4),
     &         nexc(j+5),nexc(j+6),nexc(j+7),nexc(j+8),nexc(j+9)
        enddo
C
        i2=1
        do i=1,natom
          i1=nexc(i)
          nexc(i)=i2
          i2=i2+i1
        enddo
        nexc(natom+1)=i2
C non-bonded parm
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=(ntyp**2-1)/10+1
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) nbparm(j),nbparm(j+1),nbparm(j+2),nbparm(j+3),
     &         nbparm(j+4),nbparm(j+5),nbparm(j+6),nbparm(j+7),
     &         nbparm(j+8),nbparm(j+9)
        enddo
C  residue pointers
        inext=5+(nres-1)/20
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nres-1)/10
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) nrpnt(j),nrpnt(j+1),nrpnt(j+2),nrpnt(j+3),
     &         nrpnt(j+4),nrpnt(j+5),nrpnt(j+6),nrpnt(j+7),
     &         nrpnt(j+8),nrpnt(j+9)
        enddo
        nrpnt(nres+1)=natom+1
C  bond force
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numbnd-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) kbnd(j),kbnd(j+1),kbnd(j+2),kbnd(j+3),kbnd(j+4)
        enddo
C  bond distance
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numbnd-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) rbnd(j),rbnd(j+1),rbnd(j+2),rbnd(j+3),rbnd(j+4)
        enddo
C  angle force
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numang-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) kang(j),kang(j+1),kang(j+2),kang(j+3),kang(j+4)
        enddo
C  angle values
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(numang-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) tang(j),tang(j+1),tang(j+2),tang(j+3),tang(j+4)
        enddo
C  dihedral force
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) kdih(j),kdih(j+1),kdih(j+2),kdih(j+3),kdih(j+4)
        enddo
C  dihedral preriodicity
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) ndih(j),ndih(j+1),ndih(j+2),ndih(j+3),ndih(j+4)
        enddo
C  dihedral phase
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) pdih(j),pdih(j+1),pdih(j+2),pdih(j+3),pdih(j+4)
        enddo
C  SCEE scale factor
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) scee(j),scee(j+1),scee(j+2),scee(j+3),scee(j+4)
        enddo
C  SCNB scale factor
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nptra-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) scnb(j),scnb(j+1),scnb(j+2),scnb(j+3),scnb(j+4)
        enddo
C  Lennard-Jones A
        inext=5+(natyp-1)/5
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(ntyp*(ntyp+1)/2-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) alj(j),alj(j+1),alj(j+2),alj(j+3),alj(j+4)
        enddo
        do i=1,(ntyp*(ntyp+1))/2
          alj(i)=alj(i)*12.d0
        enddo
C  Lennard-Jones B
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(ntyp*(ntyp+1)/2-1)/5
        do i=1,inext
          j=(i-1)*5+1
          read(20,997) blj(j),blj(j+1),blj(j+2),blj(j+3),blj(j+4)
        enddo
        do i=1,(ntyp*(ntyp+1))/2
          blj(i)=blj(i)*6.d0
        enddo
C  Bond list w/ H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(3*nbonh-1)/10
        do i=1,inext
          j=(i-1)*10+3
          i1=mod((j+1),3)
          if (i1.eq.0) i1=3
          i2=i1+1
          if (i2.eq.4) i2=1
          i3=i2+1
          if (i3.eq.4) i3=1
          read(20,996) bhlist(i1,j/3),bhlist(i2,(j+1)/3),
     &         bhlist(i3,(j+2)/3),bhlist(i1,(j+3)/3),bhlist(i2,(j+4)/3),
     &         bhlist(i3,(j+5)/3),bhlist(i1,(j+6)/3),bhlist(i2,(j+7)/3),
     &         bhlist(i3,(j+8)/3),bhlist(i1,(j+9)/3)
        enddo
        do i=1,nbonh
          bhlist(1,i)=bhlist(1,i)/3+1
          bhlist(2,i)=bhlist(2,i)/3+1
        enddo
C  Bond list w/o H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(3*nbona-1)/10
        do i=1,inext
          j=(i-1)*10+3
          i1=mod((j+1),3)
          if (i1.eq.0) i1=3
          i2=i1+1
          if (i2.eq.4) i2=1
          i3=i2+1
          if (i3.eq.4) i3=1
          read(20,996) balist(i1,j/3),balist(i2,(j+1)/3),
     &         balist(i3,(j+2)/3),balist(i1,(j+3)/3),balist(i2,(j+4)/3),
     &         balist(i3,(j+5)/3),balist(i1,(j+6)/3),balist(i2,(j+7)/3),
     &         balist(i3,(j+8)/3),balist(i1,(j+9)/3)
        enddo
        do i=1,nbona
          balist(1,i)=balist(1,i)/3+1
          balist(2,i)=balist(2,i)/3+1
        enddo
C  Angle list w/ H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(4*ntheth-1)/10
        do i=1,inext
          j=(i-1)*10+4
          i1=mod((j+1),4)
          if (i1.eq.0) i1=4
          i2=i1+1
          if (i2.eq.5) i2=1
          i3=i2+1
          if (i3.eq.5) i3=1
          i4=i3+1
          if (i4.eq.5) i4=1
          read(20,996) ahlist(i1,j/4),ahlist(i2,(j+1)/4),
     &         ahlist(i3,(j+2)/4),ahlist(i4,(j+3)/4),ahlist(i1,(j+4)/4),
     &         ahlist(i2,(j+5)/4),ahlist(i3,(j+6)/4),ahlist(i4,(j+7)/4),
     &         ahlist(i1,(j+8)/4),ahlist(i2,(j+9)/4)
        enddo
        do i=1,ntheth
          ahlist(1,i)=ahlist(1,i)/3+1
          ahlist(2,i)=ahlist(2,i)/3+1
          ahlist(3,i)=ahlist(3,i)/3+1
        enddo
C  Angle list w/o H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(4*ntheta-1)/10
        do i=1,inext
          j=(i-1)*10+4
          i1=mod((j+1),4)
          if (i1.eq.0) i1=4
          i2=i1+1
          if (i2.eq.5) i2=1
          i3=i2+1
          if (i3.eq.5) i3=1
          i4=i3+1
          if (i4.eq.5) i4=1
          read(20,996) aalist(i1,j/4),aalist(i2,(j+1)/4),
     &         aalist(i3,(j+2)/4),aalist(i4,(j+3)/4),aalist(i1,(j+4)/4),
     &         aalist(i2,(j+5)/4),aalist(i3,(j+6)/4),aalist(i4,(j+7)/4),
     &         aalist(i1,(j+8)/4),aalist(i2,(j+9)/4)
        enddo
        do i=1,ntheta
          aalist(1,i)=aalist(1,i)/3+1
          aalist(2,i)=aalist(2,i)/3+1
          aalist(3,i)=aalist(3,i)/3+1
        enddo
C  Dihedral list w/ H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(5*nphih-1)/10
        do i=1,inext
          j=(i-1)*2+1
          read(20,996) dhlist(1,j),dhlist(2,j),
     &         dhlist(3,j),dhlist(4,j),dhlist(5,j),
     &         dhlist(1,j+1),dhlist(2,j+1),dhlist(3,j+1),
     &         dhlist(4,j+1),dhlist(5,j+1)
        enddo
        do i=1,nphih
          dhlist(1,i)=dhlist(1,i)/3+1
          dhlist(2,i)=dhlist(2,i)/3+1
          if (dhlist(3,i).lt.0) then
            dhlist(3,i)=dhlist(3,i)/3-1
          else
            dhlist(3,i)=dhlist(3,i)/3+1
          endif
          if (dhlist(4,i).lt.0) then
            dhlist(4,i)=dhlist(4,i)/3-1
          else
            dhlist(4,i)=dhlist(4,i)/3+1
          endif
        enddo
        do i=2,nphih
          j=i-1
          if (iabs(dhlist(1,i)).eq.iabs(dhlist(1,j))) then
            if (iabs(dhlist(2,i)).eq.iabs(dhlist(2,j))) then
              if (iabs(dhlist(3,i)).eq.iabs(dhlist(3,j))) then
                if (iabs(dhlist(4,i)).eq.iabs(dhlist(4,j))) then
                  dhlist(2,i)=-dhlist(2,i)
                endif
              endif
            endif
          endif
        enddo
C  Dihedral list w/o H
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(5*nphia-1)/10
        do i=1,inext
          j=(i-1)*2+1
          read(20,996) dalist(1,j),dalist(2,j),
     &         dalist(3,j),dalist(4,j),dalist(5,j),
     &         dalist(1,j+1),dalist(2,j+1),dalist(3,j+1),
     &         dalist(4,j+1),dalist(5,j+1)
        enddo
        do i=1,nphia
          dalist(1,i)=dalist(1,i)/3+1
          dalist(2,i)=dalist(2,i)/3+1
          if (dalist(3,i).lt.0) then
            dalist(3,i)=dalist(3,i)/3-1
          else
            dalist(3,i)=dalist(3,i)/3+1
          endif
          if (dalist(4,i).lt.0) then
            dalist(4,i)=dalist(4,i)/3-1
          else
            dalist(4,i)=dalist(4,i)/3+1
          endif
        enddo
        do i=2,nphia
          j=i-1
          if (iabs(dalist(1,i)).eq.iabs(dalist(1,j))) then
            if (iabs(dalist(2,i)).eq.iabs(dalist(2,j))) then
              if (iabs(dalist(3,i)).eq.iabs(dalist(3,j))) then
                if (iabs(dalist(4,i)).eq.iabs(dalist(4,j))) then
                  dalist(2,i)=-dalist(2,i)
                endif
              endif
            endif
          endif
        enddo
C  Excluded atom list
        inext=2
        do i=1,inext
          read(20,*)
        enddo
        inext=1+(nnb-1)/10
        do i=1,inext
          j=(i-1)*10+1
          read(20,996) excl(j),excl(j+1),excl(j+2),excl(j+3),excl(j+4),
     &         excl(j+5),excl(j+6),excl(j+7),excl(j+8),excl(j+9)
        enddo
C
999     format(8(i8))
998     format(2(i8),3(8X),4(i8))
997     format(5(e16.8))
996     format(10(i8))
        close(20)
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

