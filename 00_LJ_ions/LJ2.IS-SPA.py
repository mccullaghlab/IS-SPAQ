import numpy as np

d0=-0.21

eps=2.3473
rho=0.0074

q0=1.0

pdat=np.loadtxt("sample.qp.gr1")
mdat=np.loadtxt("sample.qm.gr1")

dx=0.1
rmax=15.
dbins=np.arange(dx,rmax+dx,dx)
nbins=len(dbins)

rC=abs(d0)
rH=abs(1.1+d0)
rCl=np.sqrt(1.758**2+d0**2+2.*d0*1.758*np.cos(1.87937134))
if abs(d0) < 1.e-8:
    cC=1.
else:
    cC=d0/rC
cH=(1.1+d0)/rH
if abs(d0) < 1.e-8:
    cCl=np.cos(1.87937134)
else:
    cCl=(d0**2+rCl**2-1.758**2)/2./d0/rCl
qH=0.2659
qCl=0.0396*3.
qC=-qH-qCl

dipo=qH*rH*cH+qCl*rCl*cCl+qC*rC*cC
quad=qH*rH**2*(1.5*cH**2-0.5)+qCl*rCl**2*(1.5*cCl**2-0.5)+qC*rC**2*(1.5*cC**2-0.5)
octu=qH*rH**3*(2.5*cH**2-1.5)*cH+qCl*rCl**3*(2.5*cCl**2-1.5)*cCl+qC*rC**3*(2.5*cC**2-1.5)*cC

r0=np.append(pdat[:,0],[0.,0.])
gp=np.append(pdat[:,1],[1.,1.])
gp=np.log(gp,out=-1.E15*np.ones_like(gp),where=gp!=0.)
fp=np.append(pdat[:,2],[0.,0.])
fp[np.where(np.isnan(fp))]=0.
pp=np.append(pdat[:,4],[0.,0.])
pp[np.where(np.isnan(pp))]=0.

gm=np.append(mdat[:,1],[1.,1.])
gm=np.log(gm,out=-1.E15*np.ones_like(gm),where=gm!=0.)
fm=np.append(mdat[:,2],[0.,0.])
fm[np.where(np.isnan(fm))]=0.
pm=np.append(mdat[:,4],[0.,0.])
pm[np.where(np.isnan(pm))]=0.

Ep=pp/(1.-pp**2)*(3.-(6.*pp**2+pp**4-2.*pp**6)/5.)
Em=pm/(1.-pm**2)*(3.-(6.*pm**2+pm**4-2.*pm**6)/5.)

dintg=0.1
xx,yy=np.meshgrid(np.arange(-99.95,100,dintg),np.arange(0.05,25.,dintg))

nmax=len(r0)-2

fpC=np.zeros(nbins)
fmC=np.zeros(nbins)
fpLJ=np.zeros(nbins)
fmLJ=np.zeros(nbins)
upC=np.zeros(nbins)
umC=np.zeros(nbins)
upLJ=np.zeros(nbins)
umLJ=np.zeros(nbins)
for iR,R in enumerate(dbins):
    hR=R/2.
   
    rp=np.sqrt((xx+hR)**2+yy**2)
    drp=rp/dx-0.5
    irp=np.ndarray.astype(np.floor(drp),np.int)
    drp-=irp
    irp[np.where(irp>nmax-2)]=nmax
    irp[np.where(irp<0)]=0
    
    rm=np.sqrt((xx-hR)**2+yy**2)
    drm=rm/dr-0.5
    irm=np.ndarray.astype(np.floor(drm),np.int)
    drm-=irm
    irm[np.where(irm>nmax-2)]=nmax
    irm[np.where(irm<0)]=0
    
    g=np.exp(gp[irp]+drp*(gp[irp+1]-gp[irp])+gm[irm]+drm*(gm[irm+1]-gm[irm]))
    
    Ex=(xx+hR)*(Ep[irp]+drp*(Ep[irp+1]-Ep[irp]))/rp+(xx-hR)*(Em[irm]+drm*(Em[irm+1]-Em[irm]))/rm
    Ey=    yy *(Ep[irp]+drp*(Ep[irp+1]-Ep[irp]))/rp+    yy *(Em[irm]+drm*(Em[irm+1]-Em[irm]))/rm
    
    Ex+=np.divide(-(xx+hR)*(1.-1./eps)*3./4./np.pi/rho/dipo,rp**3,out=np.zeros_like(Ex),where=(irp==nmax) & (irm!=nmax))
    Ex+=np.divide( (xx-hR)*(1.-1./eps)*3./4./np.pi/rho/dipo,rm**3,out=np.zeros_like(Ex),where=(irm==nmax) & (irp!=nmax))
    Ey+=np.divide(-yy*(1.-1./eps)*3./4./np.pi/rho/dipo,rp**3,out=np.zeros_like(Ex),where=(irp==nmax) & (irm!=nmax))
    Ey+=np.divide( yy*(1.-1./eps)*3./4./np.pi/rho/dipo,rm**3,out=np.zeros_like(Ex),where=(irm==nmax) & (irp!=nmax))

    Etot=np.sqrt(Ex**2+Ey**2)
    px=np.divide(Ex,Etot,out=np.zeros_like(Etot),where=Etot!=0.)
    py=np.divide(Ey,Etot,out=np.zeros_like(Etot),where=Etot!=0.)

    Rzp=(px*(xx+hR)+py*yy)/rp
    Rzm=(px*(xx-hR)+py*yy)/rm
    
    cothE=np.divide(1.,np.tanh(Etot),out=np.zeros_like(Etot),where=Etot!=0.)
    Einv=np.divide(1.,Etot,out=np.zeros_like(Etot),where=Etot!=0.)

    c1=cothE-Einv
    c2=1.-2.*Einv*c1
    c3=cothE-3.*Einv*c2
    
    fpx=g*((fp[irp]+drp*(fp[irp+1]-fp[irp]))*(xx+hR)/rp)*yy
    fmx=g*((fm[irm]+drm*(fm[irm+1]-fm[irm]))*(xx-hR)/rm)*yy
    
    fpLJ[iR]= 2.*np.pi*rho*dintg**2*np.sum(fpx)
    fmLJ[iR]=-2.*np.pi*rho*dintg**2*np.sum(fmx)


    #dipole
    fpx= dipo*c1/rp**3*(3.*Rzp*(xx+hR)/rp-px)
    fmx=-dipo*c1/rm**3*(3.*Rzm*(xx-hR)/rm-px)
    #quadrupole
    fpx+= quad*(1.5*c2-0.5)/rp**4*((7.5*Rzp**2-1.5)*(xx+hR)/rp-3.*Rzp*px)
    fmx+=-quad*(1.5*c2-0.5)/rm**4*((7.5*Rzm**2-1.5)*(xx-hR)/rm-3.*Rzm*px)
    #octupole
    fpx+= octu*(2.5*c3-1.5*c1)/rp**5*((17.5*Rzp**3-7.5*Rzp)*(xx+hR)/rp-(7.5*Rzp**2-1.5)*px)
    fmx+=-octu*(2.5*c3-1.5*c1)/rm**5*((17.5*Rzm**3-7.5*Rzm)*(xx-hR)/rm-(7.5*Rzm**2-1.5)*px)
    #scale
    fpx*=g*yy
    fmx*=g*yy
    
    fpC[iR]= 2.*332.*q0*np.pi*rho*dintg**2*np.sum(fpx)
    fmC[iR]=-2.*332.*q0*np.pi*rho*dintg**2*np.sum(fmx)

hdx=dx*0.5
for i in range(nbins-1):
    upC[i+1]=upC[i]-hdx*(fpC[i]+fpC[i+1])
    umC[i+1]=umC[i]-hdx*(fmC[i]+fmC[i+1])
    upLJ[i+1]=upLJ[i]-hdx*(fpLJ[i]+fpLJ[i+1])
    umLJ[i+1]=umLJ[i]-hdx*(fmLJ[i]+fmLJ[i+1])

upC-=upC[nbins-1]
umC-=umC[nbins-1]
upLJ-=upLJ[nbins-1]
umLJ-=umLJ[nbins-1]

fout=open("sample.IS-SPA.dat","w")
for i in range(nbins):
    fout.write("{:6.3f} {:19.12e} {:19.12e} {:19.12e} {:19.12e} {:19.12e} {:19.12e} {:19.12e} {:19.12e}\n".format(dx*(i+1),fpLJ[i],fmLJ[i],fpC[i],fmC[i],upLJ[i],umLJ[i],upC[i],umC[i]))

fout.close()
