# Project: Komodo 
AME-4593-001: Space Systems and Misson Design Undergrad Project code repository for the Mission to Titan Komodo

- The Mission_to_Titan files are not fully ooptimized, yet should deliver a completely full output of the mission; not including any of the drag graphs. 

## Mission_to_Titan.py
Main program containg majority of calculations as well as data outputs and graphs. Includes mission Delta-V calculations using ephemris data.

- Note: Internal Julian/Sidreal time not impllimented as ther is an issue with using either int() or np.floor() to perform the same task as INT() in Matlb as described on page 250 of the Curtis book. Error is that Julian date is always calcualted off by 10.0 this oculd be with the rounding down function or a typo of the equations 5.47 & 5.48; the Boulet (1999) formulation.

## Mission_to_Titan_Falcon_Heavy.py
A modified Mission_to_Titan.py file; with change values for Isp and spacecraft mass assuming that for Today's open market the project could only afford/obtain the Felcon Heavy Rocket.
## Lamberts_velocity_Calcualtor_1.0.py
The calculator used for the ephemeris calculations.
## Saturn_Atmo_Drag_Delta_V_Calc.py
The simulation used for calcualting velocity change due to atmosphere drag from Saturn.
