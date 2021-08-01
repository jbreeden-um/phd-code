## High Relative Degree Control Barrier Functions Under Input Constraints

This paper was accepted for publication at the 2021 IEEE Conference on Decision and Control. 

The code in these two subfolders contain everything necessary to reproduce the results in the paper. That said, note that a more versatile and faster running version of much of this code can be found in the [Automatica](../Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories) folder. 

### Sphere

In [aa_RunSim.m](Sphere/aa_RunSim.m), select between case 1, 2, or 3, which correspond to the order of the cases in the paper. Some of these cases will seem to freeze as the agent gets close to the exclusion zone. This is normal, and is because the simulator is using a smaller time step in order to avoid chattering control inputs. The potential for chatter is one of the motivations for [future work](../Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories).

### Asteroid

Simply run [aa_RunSim.m](Asteroid/aa_RunSim.m) to view this experiment. Note that on your first run, the script will generate a file that's used to compute the "nominal path" around the asteroid. This file is too big to share via git, so it is saved during the first run. Once the simulation is running, rotate the figure that appears until you can see the trajectory.