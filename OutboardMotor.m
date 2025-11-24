%% Outboard Motor Calculations
clear
clc

% F115 Specs
PSP = 115;     % propeller shaft Power [hp]
rotMax = 5500; % maximum shaft rotation speed [rpm]

% SX190 Specs
draft0 = 1.25; % boat draft at wt [ft]
beam = 8.167;  % boat beam [ft]
LWL = 18; % length at waterline [ft]
spdmax = 37;   % max boat speed with max weight capacity [mph]
wt = 2370;     % boat weight [lb]
wtcap = 1600;  % weight capacity [lb]
wtmax = wt + wtcap; % total maximum weight [lb]
accelMax = 5; % maximum acceleration of boat [ft/s^2]

% Physical constants
g = 32.18;       % gravity [ft/s^2]
rho = 62.38;     % density of water [slug/ft^3]
nu = 2.03*10^-5; % kinematic viscosity of water [slug*s/ft^2]

%% Calculate forces on shaft
% Create vectors of Weight, Speed, and Acceleration values
spdmax = spdmax*(3600/5280); % convert [mph] to [ft/s]
Speed = linspace(0,spdmax); % create vector of Speed values
Weight = linspace(wt,wtmax); % create vector of Weight values
Accel = linspace(0,accelMax); % create vector of Acceleration values

wpArea = beam * LWL * 0.75; % cross section area at waterline (C_w = 0.75 assumed) [ft^2]
draft = @(Weight) Weight/(rho * g * wpArea); % height of boat below the water [ft]
wetArea = @(Weight) beam * draft(Weight) * 2; % wetted area [ft^2]

% Force calculations constants
Re = @(Speed, Weight) Speed .* wetArea(Weight)/nu; % Reynold's number [-]
C_f = @(Speed, Weight) 0.075./(log10(Re(Speed, Weight)) - 2).^2; % drag coefficient [-]
C_d = @(Speed, Weight) C_f(Speed, Weight) * 1.15;

% Force required to propel boat
Fd = @(Speed, Weight) (1/2)*rho * Speed.^2 .* C_d(Speed, Weight) .* wetArea(Weight); % force to overcome drag [lb]
Fp = @(Weight, Accel) (Weight/g) .* Accel; % force of acceleration [lb]

% Axial stress: sigma = F/A
r = 2; % radius of propeller shaft [in]
A = pi*r^2; % area of cross section of prop shaft [in^2]
sigmaD = @(Weight,Speed) Fd(Weight,Speed)/A; % maximum axial stress on prop shaft [lb/in^2]
sigmaA = @(Weight,Accel) Fp(Weight,Accel)/A; % maximum axial stress on prop shaft [lb/in^2]

% calculate shear stress P = T*omega
pwr = PSP*33000; % convert [hp] to [ft*lb/min]
torqueMax = pwr/rotMax; % torque on prop shaft [lb*in]

r = 2; % radius of propeller shaft [in]
J = pi*r^4/2; % polar moment of inertia [in^4]

% tau = (T*r)/J
tauMax = torqueMax*r/J; % maximum shear stress on prop shaft [lb/in^2]

%% Generate Approximate solutions

% compute the approximate solution
y0 = 1;
yapproxD = forwardEulerHR(sigmaD, Speed, y0);
yapproxA = forwardEulerHR(sigmaA, Accel, y0);
sigmaTot = yapproxD + yapproxA;
sigmaMax = max(sigmaTot);

% plot stress v speed
figure(1)
plot(Speed, yapproxD, 'r')
xlim([0, spdmax]);
xlabel('Speed [ft/s]'); ylabel('Axial Stress [psi]')
title('Axial Stress on Prop Shaft over Speed')
grid on

% plot stress v acceleration
figure(2)
plot(Accel, yapproxA, 'r')
xlim([0, accelMax]);
xlabel('Acceleration [ft/s^2]'); ylabel('Axial Stress [psi]')
title('Axial Stress on Prop Shaft over Acceleration')
grid on

%% 3-D plot of Force = f(Weight,Speed)
% generate a three dimensional plot of the function
figure(3)
[S, W] = meshgrid(Speed, Weight);
P = sigmaD(W, S);
surf(S, W, P, 'EdgeColor', 'interp');
xlabel('Speed [ft/s]'); ylabel('Weight [lb]'); zlabel('Stress on Shaft [psi]')
xlim([0, spdmax]); ylim([wt, wtmax]); 
title('Axial Shaft Stress as a Function of Weight and Speed')
shading interp

figure(4)
[A, W] = meshgrid(Accel, Weight);
P = sigmaA(A, W);
surf(A, W, P, 'EdgeColor', 'interp');
xlabel('Acceleration [ft/s^2]'); ylabel('Weight [lb]'); zlabel('Stress on Shaft [psi]')
xlim([0, accelMax]); ylim([wt, wtmax]); 
title('Axial Shaft Stress as a Function of Weight and Acceleration')
shading interp

%% Basic Assumption Scenarios

% Remax = spdmax*beam/nu;            % [-]
% Amax = beam * 0.5;                % [ft^2]
% C_dmax = 0.075/(log(Remax) - 2)^2; % [-]
% 
% % Force at max speed and weight
% Fdmax = (1/2)*rho * spdmax^2 * C_dmax * Amax; % force to overcome drag [lb]
% Famax = (wtmax/g) * accMax;
% Fmax = Fdmax + Famax;
% 
% % normal operating condtions
% spdnorm = 30*(3600/5280);  % normal operating speed [ft/s]
% wtnorm = wt + 120*4 + 100; % weight of boat with four people and gear [lb]
% 
% Renorm = spdnorm*beam/nu;            % [-]
% Anorm = width * draft0;              % [ft^2]
% C_dnorm = 0.075/(log(Renorm) - 2)^2; % [-]
% 
% % Force at normal speed and weight
% Fdnorm = (1/2)*rho * spdnorm^2 * C_dnorm * Anorm; % force to overcome drag [lb]
% Fanorm = (wtnorm/g) * accNorm;
% Fnorm = Fdnorm + Fanorm;

%% Find normal distributions of stress and shear
mu1 = spdmax/2; % center of distribution
mu2 = accelMax/2; % center of distribution

sigma1 = spdmax/6; % standard deviation
sigma2 = accelMax/6; % standard deviation

p1 = (1./(sigma1  *sqrt(2*pi))).*exp(-((yapproxA - mu1).^2) ./ (2*sigma1^2));
p2 = (1./(sigma2  *sqrt(2*pi))).*exp(-((yapproxD - mu2).^2) ./ (2*sigma2^2));

p1 = p1 ./ sum(p1);
p2 = p2 ./ sum(p2);

normStress1 = yapproxD .* p1;
normStress2 = yapproxA .* p2;

figure(5)
plot(Speed, normStress1, 'r')
xlabel('Speed [ft/s]'); ylabel('Axial Stress [psi]')
title('Axial Stress on Prop Shaft over Speed')
xlim([0, spdmax]);
grid on

figure(6)
plot(Accel, normStress2, 'r')
xlabel('Acceleration [ft/s^2]'); ylabel('Axial Stress [psi]')
title('Axial Stress on Prop Shaft over Acceleration')
xlim([0, accelMax]);
grid on

%% Calculate Stress Life

% characteristics of material (1035 normalized)
Sut = 72000; % ultimate tensile strength of material [psi]
Sy = 35000; % yield strength [psi]
FOS = 2; % Assuming FOS of 2

% Marin factors
SePrime = 0.5*Sut;
ka = 1;
kb = 0.879*2*r^-0.107; % because of torsion
kc = 0.59; % because of torsion
kd = 1;
ke = 0.868; % for 95% reliability
kf = 1;
Se = ka*kb*kc*kd*ke*kf*SePrime; % final fatigue strength

Sf = sigma_a / (1 - sigma_m / Sut); % finite life fatigue strength

% calculate fatigue life constants
f = 1.06 - 2.8*10^-3*Sut + 6.9*10^-6*Sut^2;
a = (f*Sut)^2/Se;
b = -(1/3)*log10(f*Sut/Se);

N = (sigmaTot/a).^(1/b);
Life = sum(N); 

%% Total Time life Calculation

dayTime = 5*60; % time used in one day = 5hr [min]
yrTime = dayTime*52; % time used per year = once per week [min]

rotTotal = yrTime .* rotMax; % total rotations in a year [-]

lifeTotal = Life ./ rotTotal; % life in years [-]
