% Author: James Parkus
% Date: August 26th, 2016
% Purpose: To calculate the burn time required for deorbiting a cubesat
% with a cold gas propulsion system

clear
% Gas Properties
% Specific Heat Ratio
% k = 1.4; % Nitrogen specific heat ratio (dimensionless)
% k = 1.41; % Hydrogen
k = 1.67; % Argon, Helium, Neon, Krypton, Xenon
% k = 1.18; % Freon 22

% Molar Masses (Kg/mol)
% Molar_Mass = 0.0140067; % Nitrogen
% mm = 0.03994; % Argon
% mm = 0.004002; % Helium
% mm = 0.0837; % Krypton
Molar_Mass = 0.1313; % Xenon
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
Exit_Angle = 15*(pi()/180); % Radians
Throat_Radius = Exit_Radius/1.4; % m
Throat_Area = pi()*Throat_Radius^2; % m^2 - Throat Area
Throat_Mach = 1;
Nozzle_Length = ((Exit_Radius-Throat_Radius)/(tan(Exit_Angle))); % meters

% Pressure Conditions
Chamber_Temperature = 273; % Kelvin - Chamber Temperature
Chamber_Pressure = 101330; % Pascals - Chamber Pressure
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
    F = (P*X+Q).^(1/P)-R*X;
    dF = (P*X+Q).^((1/P)-1)-R;
    Xnew = X - F/dF;
    diff = Xnew - X;
    X = Xnew;
end
Mach = 1/sqrt(X);

Throat_Temperature = Chamber_Temperature/(2/(k+1));

Mass_Flow = (k*Chamber_Pressure*Throat_Area)/(sqrt(k*Specific_Gas*Chamber_Temperature))*(2/(k+1))^((k+1)/(2*(k-1)));

Char_Velocity = (Throat_Pressure*Throat_Area)/Mass_Flow;

Exit_Pressure = Chamber_Pressure/(1+((k-1)/2)*Mach^2)^(k/(k-1));

Exit_Velocity = sqrt(((2*k*Specific_Gas*Chamber_Temperature)/(k-1))*(1-(1/(1+((k-1)/2)*Mach^2))));

Isp = Exit_Velocity/g0;

Mass_Ratio = exp(dV/(g0*Isp));

% Assuming that the initial mass is 2.66kg

Final_Mass = Initial_Mass/1.6678;

Propellant_Mass = Initial_Mass - Final_Mass;

System_Velocity = g0*Isp*log(Mass_Ratio);

% Cf = Exit_Velocity/Char_Velocity+(Exit_Area/Throat_Area)*(Exit_Pressure/Chamber_Pressure);

F = Mass_Flow*Exit_Velocity + Exit_Pressure*Exit_Area;

Burn = (Propellant_Mass*System_Velocity)/F;

Tank_Radius = (3*m_prop*R0*Chamber_Temperature/(4*pi()*Chamber_Pressure*Molar_Mass))^(1/3);
% Calculated_Tank_Radius = Tank_Radius;
% 
% rowinc = 1;

Assigned_Tank_Radius = 0.025;
% Assigned_Prop_Mass = 0.01;
for m_prop = 0.01:0.05:1
    Calculated_Chamber_Pressure = 3*m_prop*R0*Chamber_Temperature/(4*pi*Molar_Mass*Assigned_Tank_Radius^3);
    Optimal_Prop_Mass(rowinc,1) = Calculated_Chamber_Pressure;
    Optimal_Prop_Mass(rowinc,2) = m_prop;
    rowinc = rowinc + 1;
end
% 
% rowinc = 1;
% 
% Tank_Radius = 0.025;
% for Tank_Radius = 0.001:0.001:0.03
%     Calculated_Molar_Mass = 3*Assigned_Prop_Mass*R0*Chamber_Temperature/(4*pi()*Chamber_Pressure*Tank_Radius^3);
%     Optimal_Molar_Mass(rowinc,1) = Calculated_Molar_Mass;
%     Optimal_Molar_Mass(rowinc,2) = Tank_Radius;
%     rowinc = rowinc + 1;
% end
% 
% plot(Optimal_Prop_Mass(1:19,1),Optimal_Prop_Mass(1:19,2))
% title('Calculated Chamber Pressure vs Tank Radius')
% xlabel('Calculated Chamber Pressure')
% ylabel('Propellant Mass')

% plot(Optimal_Molar_Mass(1:30,1),Optimal_Molar_Mass(1:30,2))
% title('Calculated Molar Mass vs Tank Radius')
% xlabel('Tank Radius')
% ylabel('Molar Mass')
% end
%
% % Ae & At Incrementational Calculations
%
% rowinc = 1;
%
% for k = 0.9:0.1:1.7
%     AtdivA = M*((k+1)/2*(1+((k-1)/2)*M^2))^((k+1)/(2*(k-1))); % Equation 1
%     Eq1(rowinc,1) = k;
%     Eq1(rowinc,2) = AtdivA;
%     rowinc = rowinc + 1;
% end
%
% rowinc = 1;
%
% % Pressure & Hoop Stress Calculations
%
% for t = 0.001:0.001:0.01
%     for r = 0.01:0.01:0.06
%         P = 3*m*R*Chamber_Temperature/(4*pi()*r^3);
%         Hoop = P*r/t;
%         Stress(rowinc,1) = Hoop;
%         Stress(rowinc,2) = P;
%         Stress(rowinc,3) = r;
%         Stress(rowinc,4) = t;
%         rowinc = rowinc + 1;
%     end
% end
% %
% % figure()
% % plot(Stress(1:length(Stress),1),Stress(1:length(Stress),2))
% % title('Hoop Stress Vs Pressure')
% %
% % figure()
% % plot(Stress(1:length(Stress),1),Stress(1:length(Stress),3))
% % title('Hoop Stress Vs Radius')
% %
% % figure()
% % plot(Stress(1:length(Stress),1),Stress(1:length(Stress),4))
% % title('Hoop Stress Vs Wall Thickness')
% %
% % figure()
% % plot(Stress(1:length(Stress),2),Stress(1:length(Stress),3))
% % title('Pressure Vs Radius')
% %
% % figure()
% % plot(Stress(1:length(Stress),2),Stress(1:length(Stress),4))
% % title('Pressure Vs Wall Thickness')
% %
% % for Ae = 0.001:0.001:0.01
% %     m_dot = rho*v*Ae; % Equation 2
% %     Eq2(rowinc,1) = Ae;
% %     Eq2(rowinc,2) = m_dot;
% %     rowinc = rowinc + 1;
% % end
%
% rowinc = 1;
%
% % for At = 0.0001:0.0001:0.005
% % %     m_dot = Eq2(2,2); %% Uncomment when equation 2 is fixed
% %     Tt = (k*R*(At*Pt)^2)/(m_dot*R)^2; % Equation 3
% %     Eq3(rowinc,1) = At;
% %     Eq3(rowinc,2) = Tt;
% %     rowinc = rowinc + 1;
% % end
% %
% % rowinc = 1;
%
% for k = 0.9:0.1:1.7
%     AtdivAe = Me*((k+1)/2*(1+((k-1)/2)*Me^2))^((k+1)/(2*(k-1))); % Equation 4
%     Eq4(rowinc,1) = k;
%     Eq4(rowinc,2) = AtdivAe;
%     rowinc = rowinc + 1;
% end
%
% rowinc = 1;
%
% % for At = 0.0001:0.0001:0.005
% %     Vc = Pt*At/m_dot; % Equation 5
% %     Eq5(rowinc,1) = At;
% %     Eq5(rowinc,2) = Vc;
% %     rowinc = rowinc + 1;
% % end
%
% rowinc = 1;
% %
% % for At = 0.0001:0.0001:0.005
% %     F = At*Pe*Cf; % Equation 6
% %     Eq6(rowinc,1) = At;
% %     Eq6(rowinc,2) = F;
% %     rowinc = rowinc + 1;
% % end
%
% % for Ae = 0.001:0.001:0.01
% %     At = pi()*(sqrt(Ae/pi())-2*L*sin(alpha))^2; % Equation 7
% %     Eq7(rowinc,1) = Ae;
% %     Eq7(rowinc,2) = At;
% %     rowinc = rowinc + 1;
% % end
%
% rowinc = 1;
%
% % for k = 0.9:0.1:1.7
% %     AtdivAe = Me*(((k-1)/2)/(1+((k-1)/2)*Me^2))^((k+1)/2*(k-1)); % Equation 8
% %     Eq8(rowinc,1) = k;
% %     Eq8(rowinc,2) = AtdivAe;
% %     rowinc = rowinc + 1;
% % end
% %
% % h = length(Eq1);
% %
% % figure()
% % plot(FuelTb(1:length(FuelTb),1),FuelTb(1:length(FuelTb),2))
% % title('Burn Tie vs Propellant Mass')
% % xlabel('Burn Time')
% % ylabel('Propellant Mass')
% %
% % figure()
% % plot(Eq1(1:length(Eq1),1),Eq1(1:length(Eq1),2))
% % title('Equation 1')
% % xlabel('k')
% % ylabel('At/A')
% %
% % % figure()
% % % plot(Eq2(1:60,1),Eq2(1:60,2))
% % % title('Equation 2')
% %
% % figure()
% % plot(Eq3(1:length(Eq3),1),Eq3(1:length(Eq3),2))
% % title('Equation 3')
% % xlabel('Ae')
% % ylabel('M_dot')
% %
% % figure()
% % plot(Eq4(1:length(Eq4),1),Eq4(1:length(Eq4),2))
% % title('Equation 4')
% % xlabel('At')
% % ylabel('Tt')
% %
% % figure()
% % plot(Eq5(1:length(Eq5),1),Eq5(1:length(Eq5),2))
% % title('Equation 5')
% % xlabel('At')
% % ylabel('Vc')
% %
% % figure()
% % plot(Eq6(1:length(Eq6),1),Eq6(1:length(Eq6),2))
% % title('Equation 6')
% % xlabel('At')
% % ylabel('Force')
% %
% % figure()
% % plot(Eq7(1:length(Eq7),1),Eq7(1:length(Eq7),2))
% % title('Equation 7')
% % xlabel('Ae')
% % ylabel('At')
%
% % figure()
% % plot(Eq8(1:length(Eq8),1),Eq8(1:length(Eq8),2))
% % title('Equation 8')
% % xlabel('k')
% % ylabel('Ae/At')
%
% % %
% % % figure()
% % % plot(AeTb(1:41,2),AeTb(1:41,1))
% % % title('Burn Time Vs Exit Area')
% % % ylabel('Burn Time')
% % % xlabel('Exit Area')
% %
% % % figure()
% % % plot(AeAt(1:41,1),AeAt(1:41,2))
% % % title('Exit Area VS Throat Area')
% % % ylabel('Exit Area')
% % % xlabel('Throat Area')
% %
% % % figure()
% % % plot(AeTb(1:41,1),Thrust(1:41,1))
% % % title('Burn Time VS Thrust Force')
% % % xlabel('Burn Time')
% % % ylabel('Thrust Force')
% %
% % % figure()
% % % plot(AeAt(1:41,2),Thrust(1:41,1))
% % % title('Thrust Force Vs Ae')
% % % ylabel('Force')
% % % xlabel('Ae')
% % %
% % % figure()
% % % plot(AeAt(1:41,1),AeTb(1:41,1))
% % % title('Throat Area vs Burn Time')
% % % xlabel('Throat Area')
% % ylabel('Burn Time')