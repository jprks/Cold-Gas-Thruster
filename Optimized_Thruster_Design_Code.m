% Author: James Parkus
% Date: August 26th, 2016
% Purpose: To calculate the burn time required for deorbiting a cubesat
% with a cold gas propulsion system

clc
clear

% Specific Heat Ratio
% k = 1.4; % Nitrogen specific heat ratio (dimensionless)
k = 1.41; % Hydrogen
% k = 1.67; % Argon, Helium, Neon, Krypton, Xenon
% k = 1.18; % Freon 22
% k = 1.41; % Hydrogen

% Molar Masses (Kg/mol)
% mm = 0.0140067; % Nitrogen
% mm = 0.03994; % Argon
% mm = 0.004002; % Helium
% mm = 0.0837; % Krypton
% mm = 0.1313; % Xenon
% mm = 0.02018; % Neon
mm = 0.002016; % Hydrogen
% mm = 0.08647; % Freon 22

% Constants
R0 = 8.3144598; % KJ/KmolK - gas constant
Tc = 273; % Kelvin - Chamber Temperature
Pc = 1.01325; % bars - Chamber Pressure
m_prop = 0.3; % kg - Mass of Propellant
At = 0.0201;
dV = 122; % m/2 - Delta V
m0 = 1.3; % kg - Initial Mass
m1 = 1.0; % kg - Final Mass
g0 = 9.8066; % m/s^2 - Gravity at Sea Level
v = 122;
M = 1;
rho = 1;
alpha = 15*(pi()/180);
m = 0.019210; % Kg - Mass of Tank(Empty) 

% Calculations

R = R0/mm;

Isp = dV/(g0*log(m0/m1));

Ve = g0*Isp; %exit velocity

Me = sqrt(abs(4*k*R*Tc/((-k-1)*(Ve^2*(k-1)-2*k*R*Tc))));

rowinc = 1;

r = 0.025;
Ae = pi()*r^2;
At = Ae*Me*(((k-1)/2)/(1+((k-1)/2)*Me^2))^((k+1)/2*(k-1));
AeAt(rowinc,1) = At;
AeAt(rowinc,2) = Ae;

AtdivAe = Me*(((k-1)/2)/(1+((k-1)/2)*Me^2))^((k+1)/2*(k-1));

Pt = Pc/(2/(k+1))^(k/(k-1));

Tt = Tc/(2/(k+1));

rowinc = 1;

m_dot = Pt/(R*Tt)*sqrt(k*R*Tt)*At;
AeM(rowinc,1) = m_dot;
AeM(rowinc,2) = At;

rowinc = 1;

m_dot = AeM(rowinc,1);
Vc = (Pt*At)/m_dot;
Avcatm(rowinc,1) = Vc;
Avcatm(rowinc,2) = At;
Avcatm(rowinc,3) = m_dot;

Pe = Pc/(1+((k-1)/2)*Me^2)^(k/(k-1));

rowinc = 1;

Vc = Avcatm(1,1);
Ae = AeAt(rowinc,2);
Cf = Ve/Vc+(Ae/At)*(Pe/Pc);
CfA(rowinc,1) = Cf;
CfA(rowinc,2) = Vc;
CfA(rowinc,3) = Ae;
rowinc = rowinc + 1;

rowinc = 1;

Cf = CfA(rowinc,1);
F = At*Pc*Cf;
Thrust(rowinc,1) = F;
Thrust(rowinc,2) = Cf;
Thrust(rowinc,3) = AeAt(1,1)
rowinc = rowinc + 1;

rowinc = 1;

for m_prop = 0.3:0.1:1
    F = Thrust(1,1);
    tb = (m_prop*dV)/F;
    FuelTb(rowinc,1) = tb;
    FuelTb(rowinc,2) = m_prop;
    rowinc = rowinc + 1;
end
Rod = sqrt(Ae/pi());
Rid = sqrt(At/pi());
mn = 0.01;

L = ((Rod-Rid)/(tan(0.261799)));


% end

% Ae & At Incrementational Calculations

rowinc = 1;

for k = 0.9:0.1:1.7
    AtdivA = M*((k+1)/2*(1+((k-1)/2)*M^2))^((k+1)/(2*(k-1))); % Equation 1
    Eq1(rowinc,1) = k;
    Eq1(rowinc,2) = AtdivA;
    rowinc = rowinc + 1;
end

rowinc = 1;

% Pressure & Hoop Stress Calculations

for t = 0.001:0.001:0.01
    for r = 0.01:0.01:0.06
        P = 3*m*R*Tc/(4*pi()*r^3);
        Hoop = P*r/t;
        Stress(rowinc,1) = Hoop;
        Stress(rowinc,2) = P;
        Stress(rowinc,3) = r;
        Stress(rowinc,4) = t;
        rowinc = rowinc + 1;
    end
end

figure()
plot(Stress(1:length(Stress),1),Stress(1:length(Stress),2))
title('Hoop Stress Vs Pressure')

figure()
plot(Stress(1:length(Stress),1),Stress(1:length(Stress),3))
title('Hoop Stress Vs Radius')

figure()
plot(Stress(1:length(Stress),1),Stress(1:length(Stress),4))
title('Hoop Stress Vs Wall Thickness')

figure()
plot(Stress(1:length(Stress),2),Stress(1:length(Stress),3))
title('Pressure Vs Radius')

figure()
plot(Stress(1:length(Stress),2),Stress(1:length(Stress),4))
title('Pressure Vs Wall Thickness')
% 
% for Ae = 0.001:0.001:0.01
%     m_dot = rho*v*Ae; % Equation 2
%     Eq2(rowinc,1) = Ae;
%     Eq2(rowinc,2) = m_dot;
%     rowinc = rowinc + 1;
% end

rowinc = 1;

for At = 0.0001:0.0001:0.005
%     m_dot = Eq2(2,2); %% Uncomment when equation 2 is fixed
    Tt = (k*R*(At*Pt)^2)/(m_dot*R)^2; % Equation 3
    Eq3(rowinc,1) = At;
    Eq3(rowinc,2) = Tt;
    rowinc = rowinc + 1;
end

rowinc = 1;

for k = 0.9:0.1:1.7
    AtdivAe = Me*((k+1)/2*(1+((k-1)/2)*Me^2))^((k+1)/(2*(k-1))); % Equation 4
    Eq4(rowinc,1) = k;
    Eq4(rowinc,2) = AtdivAe;
    rowinc = rowinc + 1;
end

rowinc = 1;

for At = 0.0001:0.0001:0.005
    Vc = Pt*At/m_dot; % Equation 5
    Eq5(rowinc,1) = At;
    Eq5(rowinc,2) = Vc;
    rowinc = rowinc + 1;
end

rowinc = 1;

for At = 0.0001:0.0001:0.005
    F = At*Pe*Cf; % Equation 6
    Eq6(rowinc,1) = At;
    Eq6(rowinc,2) = F;
    rowinc = rowinc + 1;
end

for Ae = 0.001:0.001:0.01
    At = pi()*(sqrt(Ae/pi())-2*L*sin(alpha))^2; % Equation 7
    Eq7(rowinc,1) = Ae;
    Eq7(rowinc,2) = At;
    rowinc = rowinc + 1;
end

rowinc = 1;

% for k = 0.9:0.1:1.7
%     AtdivAe = Me*(((k-1)/2)/(1+((k-1)/2)*Me^2))^((k+1)/2*(k-1)); % Equation 8
%     Eq8(rowinc,1) = k;
%     Eq8(rowinc,2) = AtdivAe;
%     rowinc = rowinc + 1;
% end

h = length(Eq1);

figure()
plot(FuelTb(1:length(FuelTb),1),FuelTb(1:length(FuelTb),2))
title('Burn Tie vs Propellant Mass')
xlabel('Burn Time')
ylabel('Propellant Mass')

figure()
plot(Eq1(1:length(Eq1),1),Eq1(1:length(Eq1),2))
title('Equation 1')
xlabel('k')
ylabel('At/A')

% figure()
% plot(Eq2(1:60,1),Eq2(1:60,2))
% title('Equation 2')

figure()
plot(Eq3(1:length(Eq3),1),Eq3(1:length(Eq3),2))
title('Equation 3')
xlabel('Ae')
ylabel('M_dot')

figure()
plot(Eq4(1:length(Eq4),1),Eq4(1:length(Eq4),2))
title('Equation 4')
xlabel('At')
ylabel('Tt')

figure()
plot(Eq5(1:length(Eq5),1),Eq5(1:length(Eq5),2))
title('Equation 5')
xlabel('At')
ylabel('Vc')

figure()
plot(Eq6(1:length(Eq6),1),Eq6(1:length(Eq6),2))
title('Equation 6')
xlabel('At')
ylabel('Force')

figure()
plot(Eq7(1:length(Eq7),1),Eq7(1:length(Eq7),2))
title('Equation 7')
xlabel('Ae')
ylabel('At')

% figure()
% plot(Eq8(1:length(Eq8),1),Eq8(1:length(Eq8),2))
% title('Equation 8')
% xlabel('k')
% ylabel('Ae/At')

% %
% % figure()
% % plot(AeTb(1:41,2),AeTb(1:41,1))
% % title('Burn Time Vs Exit Area')
% % ylabel('Burn Time')
% % xlabel('Exit Area')
%
% % figure()
% % plot(AeAt(1:41,1),AeAt(1:41,2))
% % title('Exit Area VS Throat Area')
% % ylabel('Exit Area')
% % xlabel('Throat Area')
%
% % figure()
% % plot(AeTb(1:41,1),Thrust(1:41,1))
% % title('Burn Time VS Thrust Force')
% % xlabel('Burn Time')
% % ylabel('Thrust Force')
%
% % figure()
% % plot(AeAt(1:41,2),Thrust(1:41,1))
% % title('Thrust Force Vs Ae')
% % ylabel('Force')
% % xlabel('Ae')
% %
% % figure()
% % plot(AeAt(1:41,1),AeTb(1:41,1))
% % title('Throat Area vs Burn Time')
% % xlabel('Throat Area')
% ylabel('Burn Time')