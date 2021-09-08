## Simulations in MATLAB

To run these simulations, run [aa_RunSim.m](aa_RunSim.m). Feel free to change the parameters `use_cbf`, `use_lyapunov`, or `nonconvex_demo` to access the various simulation cases, or anything other parameters that strike you as interesting.

Simulation cases:
* Simple Reorientation
	* Comparison Controller: `use_cbf = 0; use_lyapunov = 1; nonconvex_demo = 0;`
	* ZohCBF Controller: `use_cbf = 1; use_lyapunov = 0; nonconvex_demo = 0;`
* Complex Reorientation
	* Comparison Controller: `use_cbf = 0; use_lyapunov = 1; nonconvex_demo = 1;`
	* ZohCBF Controller: `use_cbf = 1; use_lyapunov = 0; nonconvex_demo = 1;`
	* Combined Controller: `use_cbf = 1; use_lyapunov = 1; nonconvex_demo = 1;`
* 42 Spacecraft Simulator
	* Nominal Controller: `use_cbf = 0; use_lyapunov = 0;`
	* ZohCBF Controller: `use_cbf = 1; use_lyapunov = 0;`
	
The "42 Spacecraft Simulator" cases are in a [different folder](../42), but are listed above to clarify how the controller implemented in that folder compares to the controllers in this folder.