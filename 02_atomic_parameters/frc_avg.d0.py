#!/bin/python
import numpy as np
import MDAnalysis
import time
import math
import sys

#Timing
start = time.time()

#Input values
topo = "ADI.prmtop"
traj=[]
for i in range(30):
    traj.append("belly/ADI.belly{:02d}.xyz.nc".format(i))
rMax=12.
binSize=0.1

d0=0.09

#Read in atom type
fcrd=open("ADI.crd","r")
l=fcrd.readline()
l2=l.split()
natom=np.int(l2[0])
ntyp=np.int(l2[1])
ityp=np.zeros(natom,dtype=np.int)
l=fcrd.readline()
i=0
while l:
    l2=l.split()
    ityp[i]=np.int(l2[0])-1
    l=fcrd.readline()
    i+=1
fcrd.close()

#Read in LJ params
rminu=np.zeros(natom)
epsu=np.zeros(natom)
rminv=np.zeros(5)
epsv=np.zeros(5)
flj=open("LJparm.dat","r")
l=flj.readline()
l=flj.readline()
for i in range(natom):
    l=flj.readline()
    l2=l.split()
    rminu[i]=np.float(l2[0])
    epsu[i]=np.float(l2[1])
l=flj.readline()
for i in range(5):
    l=flj.readline()
    l2=l.split()
    rminv[i]=np.float(l2[0])
    epsv[i]=np.float(l2[1])
flj.close()

A=12.*np.sqrt(epsu[:,None]*epsv[None,:])*((rminu[:,None]+rminv[None,:])/2.)**12
B=12.*np.sqrt(epsu[:,None]*epsv[None,:])*((rminu[:,None]+rminv[None,:])/2.)**6

#Open trajectories
coord = MDAnalysis.Universe(topo,traj)
H2OCoord = coord.select_atoms("resname CL3")
ionsCoord = coord.select_atoms("not resname CL3 and not resname Cl- and not resname Na+")

#Iterate
rMax2=rMax**2
rMax2p=(rMax+d0)**2
hrMax=0.5*rMax
nbins=np.int(rMax/binSize+0.01)
nH2O=np.zeros((ntyp,nbins),dtype=np.int)
fH2O=np.zeros((ntyp,nbins))
nWat=len(H2OCoord.atoms)
ilist=np.array([0,1,2,3,4])
for ts in coord.trajectory:
    dims = coord.dimensions[:3]
    hdims = dims/2.
    
#    sys.stdout.write("Progress: {0:.2f}% Complete\r".format((float(ts.frame) / float(len(coord.trajectory))) * 100))
#    sys.stdout.flush()

    for ia,a in enumerate(ionsCoord.atoms):
        it=ityp[ia]
        rhat=H2OCoord.atoms[np.arange(1,nWat,5)].positions-a.position
        r2Wat=np.einsum("ij,ij->i",rhat,rhat)
        ir=np.where(r2Wat<rMax2p)[0]

        rhat=rhat[ir]
        phat=H2OCoord.atoms[ir*5+1].positions-H2OCoord.atoms[ir*5].positions
        phat/=np.sqrt(np.einsum("ij,ij->i",phat,phat))[:,None]
        rhat+=d0*phat
        r2Wat=np.einsum("ij,ij->i",rhat,rhat)
        ir2=np.where(r2Wat<rMax2)[0]
        ir=ir[ir2]
        rhat=rhat[ir2]
        r2Wat=r2Wat[ir2]
        rWat=np.sqrt(r2Wat)
        rhat/=rWat[:,None]
        ix=np.ndarray.astype(rWat[:]/binSize,np.int)
        fWat=np.zeros_like(r2Wat)
        for ib in ilist:
            sWat=H2OCoord.atoms[ir*5+ib].positions-a.position
            s2Wat=np.einsum("ij,ij->i",sWat,sWat)
            s6Wat=s2Wat**-3
            sWat=np.einsum("ij,ij->i",sWat,rhat)
            fWat+=sWat*s6Wat*(A[ia,ib]*s6Wat-B[ia,ib])/s2Wat
        np.add.at(nH2O[it],ix,1)
        np.add.at(fH2O[it],ix,fWat)

fout=open("sample.frc.dat","w")
for i in range(ntyp):
    for j in range(nbins):
        r0=(j+0.5)*binSize
        if nH2O[i,j] == 0:
            fout.write("{:7.3f} {:19.12e} {:9d} {:3d}\n".format(r0,0.,-1,i+1))            
        else:
            fout.write("{:7.3f} {:19.12e} {:9d} {:3d}\n".format(r0,fH2O[i,j]/nH2O[i,j],nH2O[i,j],i+1))

fout.close()

#Timing    
end = time.time()
t = end - start
print("\nTotal running time: {:.2f} sec".format(t))

