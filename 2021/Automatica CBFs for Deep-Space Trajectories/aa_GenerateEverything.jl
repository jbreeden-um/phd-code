
AutonomousExecution=true

case = :Flyby
cbf = :Constant
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimConstant.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlConstant.txt", force=true);

case = :Flyby
cbf = :Variable
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimVariable.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlVariable.txt", force=true);

case = :Flyby
cbf = :Integrated
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimIntegrated2S.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlIntegrated2S.txt", force=true);

case = :Flyby
cbf = :Orbital
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimIntegrated3S.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlIntegrated3S.txt", force=true);

case = :Flyby
cbf = :NewAlpha
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimConstantNewAlpha.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlConstantNewAlpha.txt", force=true);

case = :Flyby
cbf = :NoSwitch
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimConstantNoSwitch.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlConstantNoSwitch.txt", force=true);

case = :Eros
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimEros.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlEros.txt", force=true);
Base.Filesystem.mv("Switching.txt", "Results/SwitchingEros.txt", force=true);

case = :Landing
include("RunSim.jl");
Base.Filesystem.mv("Sim.txt", "Results/SimLanding.txt", force=true);
Base.Filesystem.mv("Control.txt", "Results/ControlLanding.txt", force=true);
