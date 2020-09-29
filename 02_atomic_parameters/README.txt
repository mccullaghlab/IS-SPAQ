All the atomic parameters needed for IS-SPA are found with the following codes:

--------------------------------------
LJ solvent forces.  Run `extract_LJparm.f` first, then `frc_avg.d0.py`.

Code:
	extract_LJparm.f
		- The prmtop is found at line 31.
		- The number of solute atoms is specified at line 33.
Output:
	LJparm.dat
		- Contains the LJ parameters for like pairs to be used in calculating the forces.

Script:
	frc_avg.d0.py
		- The prmtop is found at line 12.
		- The list of trajectories is found at lines 13-5.
		- The maximum distance measured is found at line 16.
		- The bin size is found at line 17.
		- The d-value is found at line 19.
		- The crd file is found at line 22.
		- The output file name is found at line 110.

Output:
	sample.frc.dat
		- Contains the average LJ solvent forces for each atom type.

--------------------------------------
Atomic g(r) and E-field(r).  Run `SPA.poisson.omp.rc12.sym.f`, then `SPA.poisson2.omp.rc12.sym.f`, then `SPA.pol2.omp.rc12.sym.f`.

Code:
	SPA.poisson.omp.rc12.sym.f
		- Number of cpu's to run is found at line 21.
		- The input .crd file name is found at line 26.
		- Run `grep " 1$" sample.gr3 | head` in terminal and put the g(r) value of a bin with a single count in line 33.
		- The input .gr3 distribution function is at line 46.
		- The output file name is found at line 177.
		- If you want to change the output bin size - change `hmax` in lines 5 and 200 and change `dx` at line 22.
		- If you want to change the rcut - change `hmax` in lines 5 and 200 and change `rcut2` at line 23.

Output:
	sample.SPA.dat
		- The initial guess for the atomic g(r) parameters to be used by `SPA.poisson2.omp.rc12.sym.f`.

Code:
	SPA.poisson2.omp.rc12.sym.f
		- Number of cpu's to run is found at line 27.
		- The input .crd file name is found at line 34.
		- Half the value of the g(r) of a bin with a single count id found at line 31. (See note for `SPA.poisson.omp.rc12.sym.f`, above)
		- The output of `SPA.poisson.omp.rc12.sym.f` is found at line 43.
		- The input .gr3 distribution function is at line 71.
		- The output file name is found at line 210.
		- If you want to change the output bin size - change `hmax` in lines 5 and 234 and change `dx` at line 28.
		- If you want to change the rcut - change `hmax` in lines 5 and 234 and change `rcut2` at line 29.
		- Also note that the dimensions/bin sizes of the .gr3 file are hard-coded in lines 118-20.  Those will need to be changed if the gr3 has different values.

Output:
	sample.SPA2.dat
		- The best guess for the atomic g(r) parameters to be used for IS-SPA.

Code:
	SPA.pol2.omp.rc12.sym.f
		- Number of cpu's to run is found at line 26.
		- The input .crd file name is found at line 35.
		- The output of `SPA.poisson2.omp.rc12.sym.f` is found at line 55.
		- The input .gr3 distribution function is at line 82.
		- The output file name is found at line 253.
		- If you want to change the output bin size - change `hmax` in lines 5 and 278 and change `dx` at line 27.
		- If you want to change the rcut - change `hmax` in lines 5 and 278 and change `rcut2` at line 28.
		- The solvent properties of the density (\AA^{-3}), dipole moment (e\AA), and dielectric constant are coded in lines 30-2.
		- There is a chance this code can diverge in the fit.  If so - change (somehow?) the initial guess of the field found in line 66.

Input:
	q.dat
		- Contains a line per atom in solute with the charge, in Amber units.

Output:
	sample.pol2.dat
		- The best guess for the atomic electric field parameters to be used for IS-SPA.  It's actually the "reduced" electric field: (p E/T - p E_{CDD}/T)
