## Quadratic Programs for High Relative Degree Spatial Constraints and Spatiotemporal Specifications with Spacecraft Applications

This [paper](https://ieeexplore.ieee.org/document/9304162) was presented at the 2020 IEEE Conference on Decision and Control.

See a video of the resultant trajectory at [https://youtu.be/9VmAR6mQoyc](https://youtu.be/9VmAR6mQoyc).

## Post Submission/Application Notes

1. Note that the way this paper considers multiple constraints potentially leads to problems with existence and uniquenss of system trajectories. This is expressed in the code as sensitivity to certain changes in the simulation parameters. A provably correct way of handling multiple constraints is found in the [2021/Automatica](../../2021/Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories) folder. The approach for handling multiple constraints `H` used in [2021/Automatica](../../2021/Automatica%20Robust%20CBFs%20for%20Satellite%20Trajectories) will also work for the constraints `h` introduced here. Alternatively, one can use a [nonsmooth simulator](https://ieeexplore.ieee.org/document/7937882). The authors encourage users to read all literature carefully, as these existence and uniqueness issues surrounding multiple constraints arise in other papers as well.