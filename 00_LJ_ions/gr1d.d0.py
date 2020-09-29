#!/bin/python
import numpy as np
import MDAnalysis as mda
import time
import sys

start = time.time()

prmtop="sample.prmtop"
traj=["sample.run1.xyz.nc","sample.run2.xyz.nc"]

selLJ=u.select_atoms("bynum 1")
selWat=u.select_atoms("resname CL3")

AH=12.*5.49048940E+06
BH=6.*1.03578910E+03
AC=12.*8.06994610E+07
BC=6.*6.45179460E+03
ACl=12.*1.70300433E+08
BCl=6.*1.23046600E+04

d0=-0.21
rMax=25.
binSize=0.1
rho=0.0074

rMax2=rMax**2
rMax2p=(rMax+abs(d0))**2
binCount=np.int(rMax/binSize+0.001)
nWat=selWat.n_atoms

uH2O=np.zeros(binCount)

nH2O=np.zeros(binCount,dtype=np.int)
pH2O=np.zeros(binCount)
fLJH2O=np.zeros(binCount)
fCH2O=np.zeros(binCount)
ngo=0
u=mda.Universe(prmtop,traj)
for ts in u.trajectory:
    ngo+=1

    vC=selWat.atoms[np.arange(1,nWat,5)].positions-selLJ.atoms[0].position
    rC2=np.einsum("ij,ij->i",vC,vC)
    ir=np.where(rC2<rMax2p)[0]
    vC=vC[ir]
    pnow=selWat.atoms[5*ir+1].positions-selWat.atoms[5*ir].positions
    pnow/=np.sqrt(np.einsum("ij,ij->i",pnow,pnow))[:,None]
    vR=vC+d0*pnow
    rR=np.sqrt(np.einsum("ij,ij->i",vR,vR))
    ir2=np.where(rR<rMax)
    vR=vR[ir2]
    rR=rR[ir2]
    vC=vC[ir2]
    pnow=pnow[ir2]
    ir=ir[ir2]
    rC=np.sqrt(rC2[ir])

    pnow=np.einsum("ij,ij->i",pnow,vR)/rR

    vH=selWat.atoms[5*ir].positions-selLJ.atoms[0].position
    rH=np.sqrt(np.einsum("ij,ij->i",vH,vH))
    vCl1=selWat.atoms[5*ir+2].positions-selLJ.atoms[0].position
    rCl1=np.sqrt(np.einsum("ij,ij->i",vCl1,vCl1))
    vCl2=selWat.atoms[5*ir+3].positions-selLJ.atoms[0].position
    rCl2=np.sqrt(np.einsum("ij,ij->i",vCl2,vCl2))
    vCl3=selWat.atoms[5*ir+4].positions-selLJ.atoms[0].position
    rCl3=np.sqrt(np.einsum("ij,ij->i",vCl3,vCl3))

    fLJnow= np.einsum("ij,ij->i",vR,vC)*(AC/rC**6-BC)/rC**8
    fLJnow+=np.einsum("ij,ij->i",vR,vH)*(AH/rH**6-BH)/rH**8
    fLJnow+=np.einsum("ij,ij->i",vR,vCl1)*(ACl/rCl1**6-BCl)/rCl1**8
    fLJnow+=np.einsum("ij,ij->i",vR,vCl2)*(ACl/rCl2**6-BCl)/rCl2**8
    fLJnow+=np.einsum("ij,ij->i",vR,vCl3)*(ACl/rCl3**6-BCl)/rCl3**8
    fLJnow/=rR

    fCnow=332.*(0.2659*vH/rH[:,None]**3-0.3847*vC/rC[:,None]**3+0.0396*vCl1/rCl1[:,None]**3+0.0396*vCl2/rCl2[:,None]**3+0.0396*vCl3/rCl3[:,None]**3)
    fCnow=np.einsum("ij,ij->i",fCnow,vR)/rR
    
    if iq == 0 or iq == 10:
        unow = 0.5*dq*332.*(0.2659/rH-0.3847/rC+0.0396/rCl1+0.0396/rCl2+0.0396/rCl3)
    else:
        unow = dq*332.*(0.2659/rH-0.3847/rC+0.0396/rCl1+0.0396/rCl2+0.0396/rCl3)

    ix=np.ndarray.astype(rR[:]/binSize,np.int)
    np.add.at(nH2O,ix,1)
    np.add.at(uH2O,ix,unow)
    np.add.at(fLJH2O,ix,fLJnow)
    np.add.at(fCH2O,ix,fCnow)
    np.add.at(pH2O,ix,pnow)

outFile = open("sample.gr1", 'w')

dxH2O=1./(4.*np.pi*binSize*rho*ngo)
pH2O/=nH2O
fLJH2O/=nH2O
fCH2O/=nH2O
for i in range(binCount):
    r0=(i+0.5)*binSize
    outFile.write("{:7.3f} {:19.12e} {:19.12e} {:19.12e} {:19.12e} {:8d}\n".format(r0,nH2O[i]*dxH2O/r0**2,fLJH2O[i],fCH2O[i],pH2O[i],nH2O[i]))
outFile.close()

end = time.time()
t = end - start
print("\nTotal running time: {:.2f} sec".format(t))

