Code to run simulations

Code:
	MD.IS-SPA.multi.sym.f
		- The number of cpu's used is at line 60.
		- The simulation time step (ps) is at line 61.
		- The simulation temperature (K) is at line 62.
		- The Andersen thermostat frequency is at line 63.
		- The harmonic restraint (kcal/mol/\AA^2) is hard coded in lines 112, 118.
		- The number of MC points/atom is found at line 1786.
		- The d value is coded at line 1793.
		- The bin size of the input histograms is hard coded at lines 1214, 1284, 1787.  Also note, the `hmax` value is coded as 120 at lines 13, 1046, 1770.
		- The solvent dielectric constant and density (\AA^{-3}) are found at lines 1790-1.
		- The .crd file is found at line 1824.
		- The atomic g(r) file is found at line 1834.
		- The atomic electric field file is found at line 1846.
		- The atomic solvent LJ forces file is found at line 1856.

Input Files:
	ran3.dat
		- Random number generator seed file.
	input.dat
		- Line 1: Input prmtop.
		- Line 2: Input restart file (with velocities).
		- Line 3: ivel - if 0, then hydrogen velocities are scaled to a mass of 12.  Use 0 if the restart is from an explicit solvent simulation and 1 otherwise.
		- Line 4: Root output name.
		- Line 5: Distance (\AA) of the harmonic restraint minimum.
		- Line 6-End: Number of atoms in selection 1 of the restraint, followed by that many lines with atom indices.  Then the same for selection 2.
