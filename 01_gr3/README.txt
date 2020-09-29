This directory contains the script for calculating the distribution function around a solute pinned in space.

Script:
   gr3d.p2.d0.py
	- The d-value is hard-coded at line 190
	- The molecule is rotated into a user defined orientation in lines 216-7.  x dimension - y1 cross y2. y dimension - y1.  z - x cross y.

Inputs:
   sample.config

Output (not included):
   sample.crd
	- Defines atom types as found in the prmtop and centered molecular coordinates to be used in subsequent codes.  The types may need to changed to one's likings.
   sample.p2.gr3
        - Distribution function.  First line has the number of cells printed.  Each line has the x,y,z-coordinates, the g(r), the average polarization vector, the polarization covariance tensor (5 columns), and the total number of counts in the bin.
