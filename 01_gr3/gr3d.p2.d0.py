#!/bin/python
__author__ = 'Greg Poisson' # He wrote in I/O parts.  REX wrote the code found in the iterate subroutine.

import numpy as np
import MDAnalysis
import time
import math
import sys


# FILE VARIABLES
configFile = sys.argv[1]

psf = None
outname = None
coordDcd = None
forceDcd = None
param = None
coord = None
force = None
temperature = None
dims = None
hdims = None

debug = False
junkCounter = 0     # counter used for debugging


# DEFUALT GLOBAL VARIABLES
rMin = 0
rMax = 25
binSize = 0.1
binCount = 0
elecPermittivity = 8.854e-12      # electrical permittivity of vacuum C^2/Jm
boltzmann = 1.9872041e-3          # boltzmann constant in kcal/(K mol)
rho = 0
epsilon = []        # list of epsilon values for all non-solvent atoms
lj_rMin = []          # list of rMin values for all non-solvent atoms



# Set debug mode from config file
def setDebug(cF):
    global debug
    txt = open(cF, 'r')
    line = txt.readline()
    while line != "END CONFIG\n":
        if line == "DEBUG MODE: ON\n":
            debug = True
        line = txt.readline()

# Get name of PSF file from config file
def getPsf(cF):
    global psf
    print('\n\t***DCD Analysis***')
    if debug:
        print('\t\tDebug Mode ON')
    else:
        print('\t\tDebug Mode OFF')
    txt = open(cF, 'r')
    while psf is None:
        line = txt.readline()
        if line == 'PSF FILE:\n':
            psf = txt.readline()[:-1]
            if debug:
                print('PSF File: {}'.format(psf))
        elif line == 'END CONFIG FILE\n':
            print('No PSF file found in config.')
            break

# Get name of PSF file from config file
def getOut(cF):
    global outname
    txt = open(cF, 'r')
    while outname is None:
        line = txt.readline()
        if line == 'OUTPUT FILE NAME:\n':
            outname = txt.readline()[:-1]
            if debug:
                print('OUTPUT File: {}'.format(outname))
        elif line == 'END CONFIG FILE\n':
            print('No OUTPUT file found in config.')
            break


# Get name of Coordinate DCD files from config file
def getCoordDCDs(cF):
    global coordDcd
    coordDcd = []
    txt = open(cF, 'r')
    while len(coordDcd) == 0:
        line = txt.readline()
        if line == 'COORD DCD FILES:\n':
            line = txt.readline()
            while line != '\n':
                if line == 'END CONFIG\n':
                    print('NO DCD FILES FOUND IN CONFIG')
                coordDcd.append(line[:-1])
                line = txt.readline()
            if debug:
                print('Coordinate DCD files: {}'.format(coordDcd))

# Set coordinate max/min and binsize
def getCoordBounds(cF):
    global dims,hdims
    txt = open(cF,'r')
    line = txt.readline()
    length1 = len("MIN DISTANCE: ")
    length2 = len("MAX DISTANCE: ")
    length3 = len("BIN SIZE: ")

    global rMin, rMax, binSize, binCount

    # scan config file for coord and bin values
    while line != "END CONFIG\n":
        line = txt.readline()
        if line[:length1] == "MIN DISTANCE: ":
            rem = -1 * (len(line) - length1)
            rMin = int(line[rem:-1])
        elif line[:length2] == "MAX DISTANCE: ":
            rem = -1 * (len(line) - length2)
            rMax = int(line[rem:-1])
        elif line[:length3] == "BIN SIZE: ":
            rem = -1 * (len(line) - length3)
            binSize = float(line[rem:-1])

# Define subset of data without solvent
def parseWater():
    # select all atoms that are not water or hydrogen
    if debug:
        print("\n--Reading DCD Data--\n\t- Parsing out WATER...\n")
    global ionsCoord,ionsForce
    global H2OCoord,H2OForce
    H2OCoord = coord.select_atoms("resname CL3")
    ionsCoord = coord.select_atoms("not resname CL3 and not resname Cl- and not resname Na+")

# Initialize MD Analysis
def initMDA():
    global coord, dims, hdims, rMax, binCount, debug,force

    m = False       # Swtiched to True if user requests an rMax value greater than the system allows
    # coordinate universe
    coord = MDAnalysis.Universe(psf, coordDcd)
    dims = [coord.dimensions[0], coord.dimensions[1], coord.dimensions[2]]
    hdims = [dims[0]/2,dims[1]/2,dims[2]/2]
    rMaxLimit = np.sqrt((dims[0]**2) + (dims[1]**2) + (dims[2]**2))
    if rMax > rMaxLimit:
        rMax = rMaxLimit
        m = True

    binCount = int((rMax - rMin)/binSize)

    if debug:
        print("--Dimensions of System--")
        print("\tTotal System Dimensions: {} A x {} A x {} A".format(dims[0], dims[1], dims[2]))
        print("\tMin Interparticle Distance Considered: {} A".format(rMin))
        if m:
            print("\tMax Interparticle Distance Considered: {} A\tThe requested rMax value was bigger than the simulation box size permits, and was truncated".format(rMax))
        else:
            print("\tMax Interparticle Distance Considered: {} A".format(rMax))
        print("\tBin Size Used: {}".format(binSize))
        print("\tBin Count: {}".format(binCount))

    # Truncate solvent out of the data
    parseWater()


    # Print log data
    printLogData(debug)

# Print log data
def printLogData(d):
    if d:
        global ionsCoord
        print("--Simulation Log Info--")
        # print list of atoms being considered
        print("DCD coord universe:", len(ionsCoord), "atom(s)")
        for i in range(0, len(ionsCoord)):
            print("\t", ionsCoord[i])
        # some general log info
        print("\nNumber of time steps in coordinate trajectory:", len(coord.trajectory))

# Iterate through all pairs of particles in all simulations,
#    identifying each pair of particles, performing computations,
#    and storing the results in a data set
def iterate():
    global plots,dims,hdims
    ngo=0
    nWat=len(H2OCoord.atoms)
    d0=0.09
    hrMax=0.5*rMax
    hrMaxp=0.5*rMax+d0
    nH2O=np.zeros((binCount,binCount,binCount),dtype=np.int)
    pH2O=np.zeros((binCount,binCount,binCount,3),dtype=np.float)
    p2H2O=np.zeros((binCount,binCount,binCount,3,3),dtype=np.float)
    if debug:
        print("-- Iterating through all particle pairs in first time step to establish pair types")
    for ts in coord.trajectory:                 # Iterate through all time steps
        dims = [coord.dimensions[0], coord.dimensions[1], coord.dimensions[2]]
        hdims = [dims[0]/2.,dims[1]/2.,dims[2]/2.]
        ngo+=1

#        sys.stdout.write("Progress: {0:.2f}% Complete\r".format((float(ts.frame) / float(len(coord.trajectory))) * 100))
#        sys.stdout.flush()

# Compute radial vectors
        if ts.frame <= 1:
            rCen=np.zeros(3)
            rot=np.zeros((3,3))
            mtot=0
            for a in ionsCoord:
                rCen+=a.position*a.mass
                mtot+=a.mass
            rCen=rCen/mtot
            ionsCoord.positions-=rCen
            y1=ionsCoord.atoms[14].position-ionsCoord.atoms[6].position
            y2=ionsCoord.atoms[8].position-ionsCoord.atoms[6].position
            if np.dot(y1,y2) < 0.:
                y2*=-1
            rot[0]=np.cross(y1,y2)
            rot[0]/=np.linalg.norm(rot[0])
            rot[1]=y1
            rot[1]/=np.linalg.norm(rot[1])
            rot[2]=np.cross(rot[0],rot[1])
            print(rot)
            ionsCoord.positions=np.dot(ionsCoord.positions,rot.T)
            if debug:
                print("-- Printing the .crd file")
            outCrd = open(outname+".crd", 'w')
            ntyp=0
            i=-1
            typlist=[]
            atyp=np.zeros(len(ionsCoord),dtype=np.int)
            for a in ionsCoord:
                i+=1
                inew=1
                for jtyp in typlist:
#                    if a.type_index == jtyp[1]:
                    if a.type == jtyp[1]:
                        atyp[i]=jtyp[0]
                        inew=0
                if inew == 1:
                    ntyp+=1
#                    typlist.append([ntyp,a.type_index])
                    typlist.append([ntyp,a.type])
                    atyp[i]=ntyp
            i=0
            outCrd.write("{:4d} {:4d}\n".format(len(ionsCoord),ntyp))
            for a in ionsCoord:
#                outCrd.write("{:3d} {:12.6f} {:12.6f} {:12.6f}\n".format(atyp[i],a.position[0]-rCen[0],a.position[1]-rCen[1],a.position[2]-rCen[2]))
                outCrd.write("{:3d} {:12.6f} {:12.6f} {:12.6f}\n".format(atyp[i],a.position[0],a.position[1],a.position[2]))
                i+=1
            outCrd.close()

        # Calculate x,y,z bin
        rWat=H2OCoord.atoms[np.arange(1,nWat,5)].positions-rCen
        rWat=np.dot(rWat,rot.T)
        ir=np.where(np.all(abs(rWat)<hrMaxp,axis=1))[0]
        rWat=rWat[ir]
        pWat=H2OCoord.atoms[ir*5+1].positions-H2OCoord.atoms[ir*5].positions
	pWat/=np.sqrt(np.einsum("ij,ij->i",pWat,pWat))[:,None]
        pWat=np.dot(pWat,rot.T)
        rWat+=pWat*d0
        ir2=np.where(np.all(abs(rWat)<hrMax,axis=1))[0]
        rWat=rWat[ir2]
        pWat=pWat[ir2]
        pWat2=np.einsum('ij,ik->ijk',pWat,pWat)
        ix=np.ndarray.astype((rWat[:]+hrMax)/binSize,np.int)
        np.add.at(nH2O,tuple(ix.T),1)
        np.add.at(pH2O,tuple(ix.T),pWat)
        np.add.at(p2H2O,tuple(ix.T),pWat2)

    dxH2O=1/(binSize**3*0.0075*ngo)
    print("here")
    np.divide(pH2O,nH2O[:,:,:,None],out=pH2O,where=nH2O[:,:,:,None]!=0)
    print("there")
    np.divide(p2H2O,nH2O[:,:,:,None,None],out=p2H2O,where=nH2O[:,:,:,None,None]!=0)
    print("open file")
    outFile = open(outname+".p2.d0.09.gr3", 'w')
    outFile.write("{:12d}\n".format(binCount*binCount*binCount))
    for i in range(binCount):
        for j in range(binCount):
            for k in range(binCount):
                outFile.write("{:7.3f} {:7.3f} {:7.3f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:8d}\n".format((i+0.5)*binSize-hrMax,(j+0.5)*binSize-hrMax,(k+0.5)*binSize-hrMax,nH2O[i][j][k]*dxH2O,pH2O[i][j][k][0],pH2O[i][j][k][1],pH2O[i][j][k][2],p2H2O[i][j][k][0][0],p2H2O[i][j][k][0][1],p2H2O[i][j][k][0][2],p2H2O[i][j][k][1][1],p2H2O[i][j][k][1][2],nH2O[i][j][k]))
    outFile.close()
####

# main program
def main():

    # access global var for config file
    global configFile, pdb
    start = time.time()

    # Read config setting for debug mode
    setDebug(configFile)

    # Get name of PSF file from config file
    getPsf(configFile)
    
    # Get name of OUTPUT file from config file
    getOut(configFile)
    
    # Get names of Coord DCD files from config file
    getCoordDCDs(configFile)


    # Define coordinate min/max and bin size
    getCoordBounds(configFile)

    # Initialize MD Analysis
    initMDA()

    # Iterate over time steps, and perform MD calculations
    iterate()

    end = time.time()
    t = end - start
    print("\nTotal running time: {:.2f} sec".format(t))

# Main program code
main()
