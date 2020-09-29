The code that analyzes the LJ ion trajectories and integrates the equations of IS-SPA.

SCRIPT:
	gr1d.d0.py
		- Prmtop file at line 9.
		- Trajectory file(s) at line 10.
		- LJ parameters A and B for solute--solvent at lines 15-20.  These are for the force, not the energy.
		- d value (\AA) at line 22.
		- Histogram maximum and bin size at lines 23-4.
		- Solvent density at line 25.
		- Output file name at line 92.

OUTPUT:
	sample.gr1
		- Column 1: Distance (\AA)
		- Column 2: g(r)
		- Column 3: Average LJ solvent force
		- Column 4: Average C solvent force
		- Column 5: Average rhat\cdot phat
		- Column 6: Counts

SCRIPT:
	LJ2.IS-SPA.py
		- d value at line 3.
		- Solvent dielectric constant and density (\AA^{-3}) at line 5-6.
		- Charge of ions (e) at line 8.
		- Input gr1 files for the two particles at lines 10-11.
		- List of distances to compute data at line 15.
		- Domain of integration found at line 57.
		- Output file name at line 145.

		
