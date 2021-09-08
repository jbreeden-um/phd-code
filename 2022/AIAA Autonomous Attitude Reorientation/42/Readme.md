## Simulations in "42"

"42" is a spacecraft simulator developed by Eric Stoneking and published by NASA. It is primarily intended for attitude control simulations, but has many other functions. It is used in this project to validate the results of our own MATLAB simulations in an outside simulation environment.

### How to Run 42

In order to avoid copying a 300+ MB simulator, this folder only contains the files that were changed in order to produce our simulations. The simulations are based off commit `46fcd3430593e4f25a76ef6a366646530ef11ece`. You will need to download that commit or a more recent one from [the source](https://github.com/ericstoneking/42), and then modify a few of the files to incorporate our code.

Some familiarity with C will be necessary to run our simulation as is. You will have to first setup "42" and get the default simulation working before adding our modifications. Windows users, see the _Docs/install-msys.txt_ file in [the source](https://github.com/ericstoneking/42/blob/master/Docs/Install-msys.txt). You will probably encounter some bugs along the way. You will also need to get a copy of the library OSQP (at time of writing, the author is not aware of any active precompiled libraries, so you'll have to compile it yourself). 

If you have trouble or think we omitted a necessary file, contact the authors for assistance. For users who are not necessarily interested in reproducing our exact simulations but still want to understand our controller, you should be able to follow along the files in this folder without going through the simulator setup. The crucial function is `PrototypeFSW` in [42fsw.c](Source/42fsw.c). Note that the simulation date in [Inp_Sim.txt](InOut/Inp_Sim.txt) and the orbit parameters in [Orb_LEO.txt](https://github.com/ericstoneking/42/blob/master/InOut/Orb_LEO.txt) were not modified at all.

### 42 Flow

See the _Docs/_ folder in [the source](https://github.com/ericstoneking/42/blob/master/Docs/) for more details of "42's" operations. In short, you compile "42" with your intended flight controller in the [42fsw.c](Source/42fsw.c) file. All the other files make up the simulation engine. You then run the "42" executable (there are input arguments available, but we don't need any of them). The executable first reads in [Inp_Sim.txt](InOut/Inp_Sim.txt) to determine the simulation environment, then reads in the [spacecraft file](InOut/aa_SC.txt) and [orbit file](https://github.com/ericstoneking/42/blob/master/InOut/Orb_LEO.txt), and assorted other parameter files. Once the simulation starts, the executable reads in the list of [spacecraft commands](InOut/aa_Cmd.txt), which tell the controller where to point. The simulation finishes after a pre-specified duration.