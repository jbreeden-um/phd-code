## Simulations in MATLAB

To run these simulations, run [aa_RunSim.m](aa_RunSim.m). Feel free to change the parameters `control_law`, `use_cbf`, `nonconvex_demo`, or `zoh_violation_demo` to access the various simulation cases, or any other parameters that strike you as interesting.

Simulation cases:
* Simple Reorientation
	* ZohCBF Controller: `control_law = 1; use_cbf = 1; nonconvex_demo = 0;`
	* Logarithmic Barrier Controller: `control_law = 2; use_cbf = 0; nonconvex_demo = 0;`
	* Sliding Mode Controller: `control_law = 3; use_cbf = 0; nonconvex_demo = 0;`
	* Model Predictive Controller: `control_law = 4; use_cbf = 0; nonconvex_demo = 0;`
		* *Note that this controller will run without additional files from another author's repository. These files are not presently available online. See comments in `mpc_solver.m`*
* Complex Reorientation
	* Logarithmic Barrier Controller: `control_law = 2; use_cbf = 0; nonconvex_demo = 1;`
	* ZohCBF Controller: `control_law = 1; use_cbf = 1; nonconvex_demo = 1;`
	* Combined Controller: `control_law = 2; use_cbf = 1; nonconvex_demo = 1;`
* 42 Spacecraft Simulator
	* PD Controller: `control_law = 1; use_cbf = 0;`
	* ZohCBF Controller: `control_law = 1; use_cbf = 1;`
	
The "42 Spacecraft Simulator" cases are in a [different folder](../42), but are listed above to clarify how the controller implemented in that folder compares to the controllers in this folder.