% This is the first part of the simulation, to determine the strength of the magents required. 
% Based on "Passive Magnetic Attitude Control for CubeSat Spacecraft"
% This section is purely for intro calculations. 
% Use this function for exporting to csv: https://www.mathworks.com/help/releases/R2019b/matlab/ref/writematrix.html

% to interpret this data, 
% use m_bar as the recommended strength for the main magnet. 
% 


close ALL
clear 

%% General Torques acting on the cubesat
% values are from the PDF documents, subject to change depending on altitude. Can be variable. 
T_aero = 8E-8; % aerodynamic
T_grav = 6E-8; % gravity gradient
T_radi = 1E-8; % radiometric torques

% to get the RMS values, for T. 
T_rms = sqrt((1/3)*(T_aero^2 + T_grav^2 + T_radi^2)); 


% Suggested minimum Magnet strength. 
Bmin = 2.0*10^-5; % Unit is Tesla, min field strength at 600km altitude. 
Bmax = 7*(pi/180); % accuracy, in degrees
m_bar = 10*(T_rms/(Bmin*sin(Bmax))) %recommended bar strength based on modified version of Santoni et Zelli. 


%% Suggested bar_main magnet design
% needs to do multiple entries and rows before figuring out m values.
Ixx = 3.6*(10^-3); %unit is kg*m^2  // This is an assumption, to be changed after getting a CAD model data. 
Iyy = 1.7*(10^-2); %same unit as Ixx. 
% Izz isn't included as it is assumed it is the same as Iyy due to the
% design of the cubic satellite. 
Beq = 2.3*(10^-5); %unit is Tesla, magnetic flux density at equator. 

G = 6.674*(10^-11); % grav constant. 
M = 5.972*(10^24); % mass of the earth, in Kg.  
m_earth = 3; % mass of cube satellite, Unit is Kg. 
a = 1.496*(10^11); % earth semi major axis, taken from online sources. 1AU = 1.496*(10^11)m. 
d = 86.4*(10^3); % seconds per one day
no = (d/2*pi)*(sqrt(G*(M+m_earth)/(a^3))); 

k1 = 0; % these values are for an iterative function, they represent the variable K. 
k2 = 28;
dk =  0.5; 

k = k1:dk:k2; 

n = 2.63*(k.^2) - 0.49 + 0.51*(Ixx/Iyy);
n_perp = 2.63*(k.^2) - 4.25+1.25*(Ixx/Iyy); 

m_res = (Iyy*n*(no^2))/(Beq); % magnetic moment of the 
m_res_perp = (Iyy*n_perp*(no^2))/(Beq);

M = [k ; n; n_perp; m_res; m_res_perp];

plot(k,m_res,'r',k,m_res_perp,'b')
xlabel('k') % label the x-axis
ylabel('m_res') % label the y-axis

csvwrite('magnetic_moment_results.csv',M)  % to export the data for analysis. 
saveas(gcf,'Plot_magnet_strength.png') % this will export the plot as a png file. 
%% Hysteresis Rods design 
%purpose of these rods are to be mounted in pairs, orthogonal to the main magnet to benefit from maximum dampening. 
%This section will determine how strong each rod should be, so that they can be effective. 

L = 80/1000; % length of H. rods, max length should be 9cm, as it is a phyiscal restriction, for a cube satellite. 80mm was picked as it could fit within the ACS module, and held with mounting brackets.
D = 25/1000; % diameter of rods, in meters. 

Hs = 100; % uni A/m, material saturation field strength, depends on material. Using HyMu-80 as example.
H = Hs; % temporary putting at as HS for now. 
uo = 1.25663706*(10^-6); %unit m*kg/((s^2)(A^2))permeability of free space, for sake of this sim, this value is considered as a constant. 
 
Vol_hyst = L*(0.25*pi*(D^2));  % volume of the hysteresis rods. 

N = ((L/D) * (4/sqrt(pi) ) +2)^(-1); % demagnetizing factor of the hyst. rods. 

U_hyst = 1.5*(10^4); % unit is H/m, value is varies depending on magnetic material, follows hysteresis loop diagram of B and H. 
U_hyst_appt = (U_hyst)/(1+(N*U_hyst)); % apparent relative permeability of the hyst. rod. 

B_hyst = uo*U_hyst_appt*H; % induced magnetic flux. Definition: The magnetic flux through a surface is the component of the magnetic field passing through that surface
B_hyst_sat = uo*U_hyst_appt*Hs; % induced magnetic flux. Definition: The magnetic flux through a surface is the component of the magnetic field passing through that surface

m_hyst_rod = (B_hyst * Vol_hyst) / uo;  % magentic moment for the rods. Definition: Magnetic Moment is defined as magnetic strength and orientation of a magnet or other object that produces a magnetic field. 
