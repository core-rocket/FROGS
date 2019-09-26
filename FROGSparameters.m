function [] = FROGSparameters
% FROGS
% ver1.2 (171107created)
%
% Parameter Setting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global l lcg0 lcgf lcgp lcp d m0 mf I0 If Ip0
global Cd Cnalpha Vpara1 Vpara2 Hpara
global LeleDeg LazDeg lLnchr
global WindModel dt n Cdv Zr WazDeg Vwaz

% length [m]
l        = 1.769;       % total length
lcg0     = 0.968;       % center of gravity @ take-off
lcgf     = 0.945;       % center of gravity @ engine-cut-off
lcgp     = 1.299;       % center of gravity of fuel & N2O
lcp      = 1.180;       % center of pressure
d        = 0.1227;      % outer diameter

% weight [kg]
 m0      = 6.601;       % weight @ take-off
 mf      = 6.000;       % weight @ engine-cut-off
%m0 = 6.0; mf = m0 - (7.197-6.790);

% moment of inetia [kgm^2]
I0       = 1.440;       % moment of inertia @take-off
If       = 1.288;       % moment of inertia @ engine-cut-off
Ip0      = 0.205;       % moment of inertia of fuel & N2O

% coefficient [-]
Cd       = 0.42;        % drag coefficient original 0.45 840m
Cnalpha  = 6.49;        % normal force coefficient

% parachute
Vpara1   = 8.5;         % falling velocity of 1st parachute [m/s]
Vpara2   = 8.5;         % falling velocity of 2nd parachute[m/s]
Hpara    = 300;         % 2nd parachute's deployment altitude  [m]

% launcher
LeleDeg  = 70;          % angle of elevation (vertical=90deg) [deg]
LazDeg   = 161;         % azimuth (east=0deg / south=270deg) [deg]
lLnchr   = 5.000;       % length [m]

% wind
WindModel= 1;           % model of wind speed
                        % 1:power model / 2:uniformity
Cdv      = 6.0;         % coefficient [-]
                        % in case of WindModel=1
WazDeg   = 161;         % azimuth of wind [deg]
                        % east=0deg / south=270deg
Vwaz     = 3;           % wind speed [m/s]
Zr       = 5;           % altitude anemometer located  [m]

% simulation
dt       = 0.01;        % simulation step [s]
                        % MUST BE 0.01s OR LESS!!!
n        = 10000;       % maximum nunber of simulation steps [-]
                        % if the rocket doesn't reach the ground,
                        % change 'n' to bigger one.