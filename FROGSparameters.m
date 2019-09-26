function [] = FROGSparameters
% FROGS
% ver1.4 (181124created)
%
% Parameter Setting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global l lcg0 lcgf lcgp lcp d m0 mf I0 If 
global Cd Cnalpha Vpara1 Vpara2 Hpara
global LeleDeg LazDeg lLnchr
global WindModel dt n Cdv Zr WazDeg Vwaz
global PMT

FROGSparain;        % parainの読み込み
% length [m]
l        = PMT(1);       % total length
d        = PMT(2);       % outer diameter
lcg0     = PMT(3);       % center of gravity @ take-off
lcgf     = PMT(4);       % center of gravity @ engine-cut-off
lcgp     = PMT(5);       % center of gravity of fuel & N2O

% weight [kg]
m0       = PMT(6);       % weight @ take-off
mf       = PMT(7);       % weight @ engine-cut-off
%m0 = 6.0; mf = m0 - (7.197-6.790);

% moment of inetia [kgm^2]
I0       = PMT(8);       % moment of inertia @take-off
If       = PMT(9);       % moment of inertia @ engine-cut-off

% coefficient [-]
lcp      = PMT(10);       % center of pressure
Cd       = PMT(11);        % drag coefficient original 0.45 840m
Cnalpha  = PMT(12);       % normal force coefficient

% parachute
Vpara1   = 12;           % falling velocity of 1st parachute [m/s]
Vpara2   = 8;           % falling velocity of 2nd parachute[m/s]
Hpara    = 300;           % 2nd parachute's deployment altitude  [m]

% launcher
LeleDeg  = 70;          % angle of elevation (vertical=90deg) [deg]
LazDeg   = 161;         % azimuth (east=0deg / south=270deg) [deg]
lLnchr   = 5.0;           % length [m]

% wind
WindModel= 1;           % model of wind speed
                        % 1:power model / 2:uniformity
Cdv      = 6.0;         % coefficient [-]
                        % in case of WindModel=1
WazDeg   = 0;           % azimuth of wind [deg]
                        % east=0deg / south=270deg
Vwaz     = 3;           % wind speed [m/s]
Zr       = 5;          % altitude anemometer located  [m]

% simulation
dt       = 0.01;        % simulation step [s]
                        % MUST BE 0.01s OR LESS!!!
n        = 10000;       % maximum nunber of simulation steps [-]
                        % if the rocket doesn't reach the ground,
                        % change 'n' to bigger one.