## High Relative Degree Control Barrier Functions Under Input Constraints

This paper was accepted for publication at the 2021 IEEE Conference on Decision and Control. 

The code in these two subfolders contain everything necessary to reproduce the results in the paper. That said, note that a more versatile and faster running version of much of this code can be found in the [Automatica](../Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories) folder. 

### Sphere

In [aa_RunSim.m](Sphere/aa_RunSim.m), select between case 1, 2, or 3, which correspond to the order of the cases in the paper. Some of these cases will seem to freeze as the agent gets close to the exclusion zone. This is normal, and is because the simulator is using a smaller time step in order to avoid chattering control inputs. The potential for chatter is one of the motivations for [future work](../Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories).

### Asteroid

Simply run [aa_RunSim.m](Asteroid/aa_RunSim.m) to view this experiment. Note that on your first run, the script will generate a file that's used to compute the "nominal path" around the asteroid. This file is too big to share via git, so it is saved during the first run. Once the simulation is running, rotate the figure that appears until you can see the trajectory. See a video of the resultant trajectory at [https://youtu.be/JKj3PUrYnEg](https://youtu.be/JKj3PUrYnEg).


## Post Submission/Application Notes

### Sphere

Note that this paper emphasizes what can be done, not necessarily what should be done. As such, the CBFs use an arbitrary class-K function that is not tuned to this particular problem. Because of this lack of tuning, the simulations are subject to errors from effects not considered, namely the use of numerical integration instead of a truly continuous control law. We compensated for this by using a careful integration method (see [ode11.m](ode11.m)) and small time-step. We note two somewhat obvious improvements:
1. A solver designed for stiff ODEs might produce smoother results with larger time-steps. We have had mixed results with this in subsequent work.
2. Better class-K functions could be used. The choices `alpha = -0.4*H` or `alpha = -8*H*abs(H)` allow the code to run with the current integrator and a constant setting of `RelTol=1e-2` such that the entire simulation finishes much faster.

Think of running the code as written as analogous to running a high-bandwidth controller. Unsurprisingly, a small time-step is needed. Using a less-steep class-K function is just as much provably safe, but results in the trajectory approaching the unsafe region slower, so it is not affected as greatly by numerical integration error. 

In practice, the authors recommend that users always try to tune their class-K functions to the problem at hand to avoid needing very small time-steps. Tuning also prevents chattering control inputs. However, in the spirit of truly verifying that this method is valid **for any class-K function**, we chose an arbtriary class-K function (which turned out to be "high-bandwidth"), and adjusted our numerical integration method in the code in this folder until it worked for this class-K function. Alternatively, one can directly compensate for numerical integration effects, as was partially done in our [L-CSS paper](../L-CSS%20CBFs%20for%20Sampled%20Data%20Systems).

### Asteroid

While the current asteroid simulation only considers a subset of the mesh points at any time step representing a conical sensor region (the red points in the video), the current method is not provably correct (and thus not elaborated upon in the paper). A provably correct method is presented in the [Automatica](../Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories) folder. 