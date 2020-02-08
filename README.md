# Robust-Attitude-Controller-for-UAV-Using-Dynamic-Inversion-and-Extended-State-Observer-controller
A robust feedback linearization controller is presented for attitude control of an unmanned aerial vehicle (UAV). The objective of this controller is to make the roll angle, pitch angle, and yaw angle track the given trajectories(commands) respectively. This design is developed using dynamic inversion and extended state observer (ESO). Firstly, dynamic inversion is used to linearize and decouple UAV attitude system into three single-input-single-output (SISO) systems, then three proportional-derivative (PD) controllers are designed for these linearized systems. Extended state observers are used to estimate and compensate unmodeled dynamics and extent disturbances. Simulation results show that the proposed controller is effective and robust.
 this Matlab implementation of this artcicle  https://www.researchgate.net/publication/241163325_Robust_Attitude_Controller_for_Unmanned_Aerial_Vehicle_Using_Dynamic_Inversion_and_Extended_State_Observer