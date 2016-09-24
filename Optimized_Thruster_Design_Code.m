% Author: James Parkus
% Date: August 26th, 2016
% Purpose: To calculate the burn time required for deorbiting a cubesat
% with a cold gas propulsion system
clc
clear
% Gas Properties
% Specific Heat Ratio
k = 1.4; % Nitrogen specific heat ratio (dimensionless)
% k = 1.41; % Hydrogen
% k = 1.67; % Argon, Helium, Neon, Krypton, Xenon
% k = 1.18; % Freon 22

% Molar Masses (Kg/mol)
Molar_Mass = 0.0140067; % Nitrogen
% mm = 0.03994; % Argon
% mm = 0.004002; % Helium
% mm = 0.0837; % Krypton
% Molar_Mass = 0.1313; % Xenon
% mm = 0.02018; % Neon
% mm = 0.002016; % Hydrogen
% wmm = 0.08647; % Freon 22
% Nozzle Properties
% Constants
rowinc = 1;
R0 = 8.3144598; % J/molK - gas constant
m_prop = 0.3; % kg - Mass of Propellant
dV = 122; % m/2 - Delta V Required
Initial_Mass = 2.66; % kg - Final Mass - GUESS
g0 = 9.8066; % m/s^2 - Gravity at Sea Level
v = 122;
m = 0.019210; % Kg - Mass of Tank(Empty)
Specific_Gas = R0/Molar_Mass; % J/molK -

% Nozzle Geometry & Throat Conditions
Exit_Radius = 0.025; % meters
Exit_Area = pi()*(Exit_Radius)^2; % meters ^2
Exit_Angle = 15*(pi/180); % Radians
Throat_Radius = Exit_Radius/1.4; % m
Throat_Area = pi()*Throat_Radius^2; % m^2 - Throat Area
Throat_Mach = 1;
Nozzle_Length = ((Exit_Radius-Throat_Radius)/(tan(Exit_Angle))); % meters

% Pressure Conditions
Chamber_Temperature = 273; % Kelvin - Chamber Temperature
Chamber_Pressure = 970000; % Pascals - Chamber Pressure
Throat_Pressure = Chamber_Pressure*(2/(k+1))^(k/(k-1));

% Calculations
AtVAe = Throat_Area/Exit_Area; % dimensionless

% Mach Calculation
A = Exit_Area;
At = Throat_Area;
P = 2/(k+1);
Q = 1-P;
R = (A/At).^((2*Q)/P);
a = Q.^(1/P);
r = (R-1)/(2*a);
X = 1/((1+r)+sqrt(r*(r+2)));  % Initial Guess
diff = 1;
while abs(diff) > .0001
    Thrust_Force = (P*X+Q).^(1/P)-R*X;
    dF = (P*X+Q).^((1/P)-1)-R;
    Xnew = X - Thrust_Force/dF;
    diff = Xnew - X;
    X = Xnew;
end
Mach = 1/sqrt(X);

%Newton-Raphson Approxiamtion method courtesy of Phil Lindin's Thruster
%Design Utility

Throat_Temperature = Chamber_Temperature/(2/(k+1));

Mass_Flow = (k*Chamber_Pressure*Throat_Area)/(sqrt(k*Specific_Gas*Chamber_Temperature))*(2/(k+1))^((k+1)/(2*(k-1)));

Char_Velocity = (Throat_Pressure*Throat_Area)/Mass_Flow;

Exit_Pressure = Chamber_Pressure/(1+((k-1)/2)*Mach^2)^(k/(k-1));

Exit_Velocity = sqrt(((2*k*Specific_Gas*Chamber_Temperature)/(k-1))*(1-(1/(1+((k-1)/2)*Mach^2))));

Isp = Exit_Velocity/g0;

Mass_Ratio = exp(dV/(g0*Isp));

% Assuming that the initial mass is 2.66kg

Final_Mass = Initial_Mass/Mass_Ratio;

Propellant_Mass = Initial_Mass - Final_Mass;

System_Velocity = g0*Isp*log(Mass_Ratio);

% Cf = Exit_Velocity/Char_Velocity+(Exit_Area/Throat_Area)*(Exit_Pressure/Chamber_Pressure);

Thrust_Force = Mass_Flow*Exit_Velocity + Exit_Pressure*Exit_Area;

Burn = (Propellant_Mass*System_Velocity)/Thrust_Force;

% Tank_Volume = (Molar_Mass*R0*Chamber_Temperature)/(Chamber_Pressure*Assigned_Propellant_Mass);

% Calculated_Tank_Radius = Tank_Radius;
% 
rowinc = 1;

Assigned_Tank_Radius = 0.025;
% Assigned_Propellant_Mass = 0.5; % kg

% for Chamber_Pressure = 700000:10000:7000000
Tank_Radius = (3*Molar_Mass*R0*Chamber_Temperature/(4*pi()*Chamber_Pressure*Propellant_Mass))^(1/3);
%     Optimal_Tank_Radius(rowinc,1) = Tank_Radius;
%     Optimal_Tank_Radius(rowinc,2) = Chamber_Pressure;
%     rowinc = rowinc + 1;
% end

% plot(Optimal_Tank_Radius(1:201,1),Optimal_Tank_Radius(1:201,2))
% xlabel('Tank Radius');
% ylabel('Chamber Pressure')
% Calculated_Tank_Radius = (3*Tank_Volume/(4*pi))^(1/3);

% % Assigned_Prop_Mass = 0.01;
% for m_prop = 0.01:0.05:1
%     Calculated_Chamber_Pressure = 3*Molar_Mass*R0*Chamber_Temperature/(4*pi*m_prop*Assigned_Tank_Radius^3);
%     Optimal_Prop_Mass(rowinc,1) = Calculated_Chamber_Pressure;
%     Optimal_Prop_Mass(rowinc,2) = m_prop;
%     rowinc = rowinc + 1;
% end




























