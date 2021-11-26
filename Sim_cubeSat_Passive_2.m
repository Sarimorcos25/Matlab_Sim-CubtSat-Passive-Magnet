% This is the second part of the simulation, where the focus is on the simulated perfumance of the cube-sat. 
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
wx0 = 10;
wy0 = 5;
wz0 = 5;
% intial position of satellite. 
thetax0 = 0;
thetay0 = 0;
thetaz0 = 0;

% Finding the external torque for the main bar magnet. 
Heq =  18.3; %Unit: A/m, this the magnetic field strength at 600km. (from paper, can change depending on the altitude)

u = 0:1:179;
i = 55 * (3.14/180); %this is rads, 

% these equations are the magnetic field strength, in ECEF direction vector (earth-centered, earth-fixed). 
H1 = 3*Heq*sin(i)*cos(i)*(sin(u).^2);
H2 = -3*Heq.*sin(i).*sin(u).*cos(u);
H3 = Heq*(1-3*(sin(i).^2) * (sin(u)).^2); 

plot(u,H1,'r',u,H2,'B',u,H3,'G')
xlabel('Argument of Latitude, u') % label the x-axis
ylabel('magnetic strengh in ECEF directions, H "A/m" ') % label the y-axis
saveas(gcf,'Magnetic_strength VS latitude.png') % this will export the plot as a png file.

% for the simulation of the hysteresis rods, 
p = (1/Mat_Hc)*tan( (Mat_Br*22/7)/(2*Mat_Bs));

H1 = -30; 
H2 = 30;
dH = 0.1;
H = H1:dH:H2; % for plotting purposes.  
Bhyst_1 = (14/22)*Mat_Bs*atan(p*(H + Mat_Hc) );
Bhyst_2 = (14/22)*Mat_Bs*atan(p*(H - Mat_Hc) );

Bhyst = [H; Bhyst_1; Bhyst_2]; % the equation has a +/- symbol. this is a way to write the equations in matlab. 
csvwrite('Hysteresis_response_sim.csv',Bhyst)

plot(H,Bhyst_1,'r',H,Bhyst_2,'b')
xlabel('magnetic field strength, H [A/m]') % label the x-axis
ylabel('Rod magnetic flux density Bhyst [Am2]') % label the y-axis

saveas(gcf,'Hysteresis simulation.png') % this will export the plot as a png file.
% body fixed angular acceleration. Advisable to use ODE45, ODE23 can be used 
a_x = (1/Ixx)*[Lx - (Izz-Iyy)*wy*wz];
a_b = (1/Iyy)*[Ly - (Ixx-Izz)*wz*wx];
a_c = (1/Izz)*[Lz - (Iyy-Izz)*wx*wz];
