% This is the third part of the simulation, where the focus is on simulating the expected rotation/velocity of the satellite 
% Based on "Passive Magnetic Attitude Control for CubeSat Spacecraft"

close ALL
clear 
% the values shown here are based on the magentic material 
%material chosen: HyMu-80
OD = 0.1;% cm, diameter of rod
LG = 15;% cm, length of rod.
Mat_Hc = 0.96;% A/m
Mat_Br = 0.35;% Tesla
Mat_Bs = 0.74;% Tesla, 
App_Hc = 12;% A/m
App_Br = 0.004;% Tesla
App_Bs = 0.025;% Tesla, 

% mass moment of inertia, for the satellite 
Ixx = 3.6*(10^-3); %unit is kg*m^2  // This is an assumption, can be changed if there is a 3D cad model from the whole satellite.  
Iyy = 1.7*(10^-2); %same unit as Ixx. 
Izz = Iyy; % Izz treated the same as Iyy as the design of the cube sat. showed it. 

% intial angular velocities

w_initial = [10 5 5]; % [wx wy wz]

% intial position of satellite. 
thetax0 = 0;
thetay0 = 0;
thetaz0 = 0;

a_x = (1/Ixx)*[Lx - (Izz-Iyy)*wy*wz];
a_b = (1/Iyy)*[Ly - (Ixx-Izz)*wz*wx];
a_c = (1/Izz)*[Lz - (Iyy-Izz)*wx*wz];

% there is no t, though. How do I fix this? to avoid the t? 
function Fw = Sim_cubeSat_Passive_3(t,w);
    Fw(1,1) = (1/Ixx)*[Lx - (Izz-Iyy)*wy*wz];
    Fw(2,1) = (1/Iyy)*[Ly - (Ixx-Izz)*wz*wx];
    Fw(3,1) = (1/Izz)*[Lz - (Iyy-Izz)*wx*wz];

% this is to use the ODE45 solver. 
[tw, wi] = ode45('Fw', w_initial, 1)