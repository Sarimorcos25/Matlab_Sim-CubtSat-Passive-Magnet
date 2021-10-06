% This is the second part of the simulation, where it focus on the simulated perfumance of the cube-sat. 
% Based on "Passive Magnetic Attitude Control for CubeSat Spacecraft"


close ALL
clear 

% mass moment of inertia, for the satellite 
Ixx = 3.6*(10^-3); %unit is kg*m^2  // This is an assumption, to be changed after getting a CAD model data. 
Iyy = 1.7*(10^-2); %same unit as Ixx. 
Izz = Iyy; % Izz treated the same as Iyy as the design of the cube sat. showed it. 


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

% for the simulation of the hysteresis rods, 
