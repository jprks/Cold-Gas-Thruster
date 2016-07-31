% Author: James Emerson Parkus
% Date: June 28th, 2016
% Purpose: To increment unknown variable magnitudes in conjuction with the
% known variable magnitudes to present a range of values for the design
% team to consider

clc
clear

%%% Known Variables
DeltaV = 122; %m/
Isp = 0; % seconds
gc = 9.81; %m/s^2
FuelWeight = 0.300; %kg
R = 0.008314; % KJ/KgK, MUST BE CHANGED TO ACCODMATE ALTITUDE
%R = PV/nT - Equation needs to be solved incrementally
gamma = 1.400; % KJ/KgK 
mass0 = 1.3; %kg
mass1 = 1.0; %kg
T = 263.15; % Temperature of gas, Kelvin (Increment to )
r = 0.001; % internal of the tank, m
A = 0.014975; % area m^2
Ae = 0.014975; %area m^2
t = 0.001; % meters
d = 1.251; % Kg/m^3
Vs = 343; % 343 m/s
m = 0.019210; % Kg
Vol = 0.000007; % m^3
s = 0.014975; % Exit Area, m^2 - Might want to run dimensional analysis on this
% value because the paper may have a typo

% Arrays
rowinc = 1;
datainc = 0;
datainc2 = 0;
Isp_Check = (1);
IspS = (1);
Stresst = (1);
tcheck = (1);
ThrustForce = (1);
ExitV = (1);
ExitP = (1);
ExitM = (1);
Chamber= (1);
PressureT = (1);
Sonic = (1);
CharV = (1);
AeA = (1);
% Isp Calculation

Isp = (DeltaV)/(gc*(log(mass0/mass1)));

% Isp Calculation check

rowinc = 1;

for IspLoop = 30:1:50
    V = gc*IspLoop*log(mass0/mass1);
    Isp_Check(rowinc, 1) = V;
    Isp_Check(rowinc, 2) = IspLoop;
    rowinc = rowinc + 1;
end

% Fuel Weight check

FuelWeightCheck = mass1*[exp((DeltaV)/(gc*Isp))-1];

% Pressure Calculation

P = (m*R*T)/Vol; % bars

% Internal Stress Calculation

rowinc = 1;

for t = 0.001:0.001:0.01
    stress = (P*r)/(2*t);
    Stresst(rowinc,1) = t;
    Stresst(rowinc,2) = stress;
    rowinc = rowinc + 1;
end

% Hoop Stress Calculation

hoop = (P*r)/t;

% Exit Velocity Calculation

rowinc = 1;

for Isp = 30:1:50
    Ve = Isp*gc;
    ExitV(rowinc,1) = Ve;
    ExitV(rowinc,2) = Isp;
    rowinc = rowinc + 1;
end

% Mass Flow Rate Calculation

rowinc = 1;

for mm = 1:21
    if Isp_Check(mm,1) < 122
        Isp = Isp_Check(mm,2);
    end
end

m_dot = FuelWeight/Isp;

% Thrust Force Calculation

rowinc = 1;
datainc = 1; 

for jj = (Ve-50):1:(Ve+50)
    Thrust = m_dot*jj;
    ThrustForce(rowinc,1) = Thrust;
    ThrustForce(rowinc,2) = jj;
    rowinc = rowinc + 1;
end

% Exit Pressure Calculation

AeA(1,1) = Ae;
rowinc = 1;

for ss = 1:1:10
    AeA(rowinc,1) = Ae;
    AeA(rowinc,2) = ss;
    Ae = Ae + 0.01;
    rowinc = rowinc + 1;
end

rowinc = 1;
datainc = 1;

for s = 1:1:10
    ss = AeA(datainc,1);
    Pe = Thrust/ss;
    ExitP(rowinc,1) = Pe;
    ExitP(rowinc,2) = ss;
    rowinc = rowinc + 1;
    datainc = datainc + 1;
end

% Exit Velocity Mach Number Calculation

rowinc = 1;
datainc = 1;

for i = 1:1:21
    ii = ExitV(datainc,1);
    Me = ii/Vs;
    ExitM(rowinc,1) = Me;
    ExitM(rowinc,2) = ii;
    rowinc = rowinc + 1;
    datainc = datainc + 1;
end

% Chamber Pressure in Nozzle Calculation
rowinc = 1;
datainc = 1;
datainc2 = 1;

for q = 1:1:10
    qq = ExitP(datainc,1);
    for w = 1:1:21
        ww = ExitM(datainc2,1);
        Pc = (qq)/(1+((gamma-1)/2)*(ww^2))^((gamma)/(gamma-1));
        Chamber(rowinc,1) = Pc;
        Chamber(rowinc,2) = ww;
        Chamber(rowinc,3) = qq;
        rowinc = rowinc + 1;
        datainc2 = datainc2 + 1;
    end
    datainc = datainc + 1;
    datainc2 = 1;
end

% Thrust Force for Infinite Expansion

% Thrust = At*Pc*gamma((2/(gamma-1))*(2/(gamma+1))*(1-(Pe/Pc))) + Pe*Ae;
% At IS UNKNOWN

% Pressure at Throat

rowinc = 1;
datainc = 1;

for e = 1:1:210
    ee = Chamber(datainc,1);
    Pt = ee*(1+((gamma-1)/2))^((gamma/(gamma-1)));
    PressureT(rowinc,1) = Pt;
    PressureT(rowinc,2) = ee;
    datainc = datainc + 1;
    rowinc = rowinc + 1;
end

% Sonic Velocity of Gas

rowinc = 1;
datainc = 1;

for t = 1:1:210
    tt = PressureT(datainc,1);
    a0 = ((gamma*tt)/d);
    Sonic(rowinc,1) = a0;
    Sonic(rowinc,2) = tt;
    rowinc = rowinc + 1;
    datainc = datainc + 1;
end

% Characteristic Velocity

rowinc = 1;
datainc = 1;

for y = 1:1:210
    yy = Sonic(datainc,1);
    C = yy/(gamma*(2/(gamma+1))^((gamma+1)/(2*(gamma-1))));
    CharV(rowinc,1) = C;
    CharV(rowinc,2) = yy;
    rowinc = rowinc + 1;
    datainc = datainc + 1;
end

figure()
plot(CharV(1:210,1),CharV(1:210,2))
title('Characteristic Velocity')
xlabel('Characteristic Velocity (m/s)')
ylabel('Sonic Velocity of the Gas (m/s)')

figure()
plot(ExitM(1:21,1),ExitM(1:21,2))
title('Mach Number')
xlabel('Mach Number')
ylabel('Exit Velocity (m/s)')

figure()
plot(ExitP(1:10,1),ExitP(1:10,2))
title('Exit Pressure')
xlabel('Exit Pressure (bars)')
ylabel('Exit Area (cm^2)')

figure()
plot(ThrustForce(1:101,1),ThrustForce(1:101,2))
title('Thrust Force')
xlabel('Thrust Force (Newtons)')
ylabel('Exit Velocity (m/s)')

figure()
plot(Stresst(1:10,1),Stresst(1:10,2))
title('Stress')
xlabel('Stress')
ylabel('Time (s)')

figure()
plot(ExitV(1:21,1),ExitV(1:21,2))
title('Exit Velocity')
xlabel('Exit Velocity (m/s)')
ylabel('Specific Impulse (s)')












