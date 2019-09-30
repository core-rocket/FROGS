function [Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4] = FROGSprset()
% FROGS
% ver1.8 (190807edited)
% for NSE15th
% Initial value setting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global g S Lele LeleDeg Laz LazDeg Waz WazDeg thrust
global the0 psi0 tThrust Cmq Cnalpha n SIMULATION 
global log_t log_T log_m log_I log_lcg lcp lcg0 l d dt
global log_rho log_Vw log_Vab log_Va log_alpha
global log_bet log_D log_Y log_N log_Fe log_Aeab
global log_Ae log_Ve log_Xe log_Kj log_Ka log_Veab
global log_omg log_q log_the log_psi m0 mf mp0
global Ip0 lcgp mm0 WazH WazDegH
% constant
g       = 9.80665;           % acceleration of gravity[m/s^2]
S       = pi*d*d/4;          
mp0     = m0-mf;
Lele    = -LeleDeg*pi/180;
Laz     = LazDeg*pi/180;
Waz     = WazDeg*pi/180;
tThrust = size(thrust,1)*dt;
Cmq     = -Cnalpha/2*((lcp-lcg0)/l)^2;
Ip0     = (lcgp-lcg0)^2*(mm0);  % moment of inertia of fuel & N2O
WazH     = WazDegH*pi/180;

% initial values
the0    = Lele;
psi0    = Laz;
THE 	= acos(((cos(the0)*cos(psi0)+cos(psi0)+cos(the0)-1))/2);
lambda  = [(-sin(the0)*sin(psi0))/(2*sin(THE));
           (sin(the0)*cos(psi0)+sin(the0))/(2*sin(THE));
           (cos(the0)*sin(psi0)+sin(psi0))/(2*sin(THE))];
q       = [lambda(1)*sin(THE/2);
           lambda(2)*sin(THE/2);
           lambda(3)*sin(THE/2);
           cos(THE/2)];
Xe      = [0;0;0];
Ve      = [0;0;0];
omg     = [0;0;0];

% size
if SIMULATION == 1 || 2
    log_t     = zeros(1,n);
    log_T     = zeros(1,n);
    log_m     = zeros(1,n);
    log_I     = zeros(1,n);
    log_lcg   = zeros(1,n);
    log_rho   = zeros(1,n);
    log_Vw    = zeros(3,n);
    log_Vab   = zeros(3,n);
    log_Va    = zeros(1,n);
    log_alpha = zeros(1,n);
    log_bet   = zeros(1,n);
    log_D     = zeros(1,n);
    log_Y     = zeros(1,n);
    log_N     = zeros(1,n);
    log_Fe    = zeros(3,n);
    log_Ae    = zeros(3,n);
    log_Ve    = zeros(3,n);
    log_Aeab  = zeros(1,n);
    log_Veab  = zeros(1,n);
    log_Xe    = zeros(3,n);
    log_Kj    = zeros(1,n);
    log_Ka    = zeros(1,n);
    log_omg   = zeros(3,n);
    log_q     = zeros(4,n);
    log_the   = zeros(1,n);
    log_psi   = zeros(1,n);
end

% output
Ve1  = Ve(1);
Ve2  = Ve(2);
Ve3  = Ve(3);
Xe1  = Xe(1);
Xe2  = Xe(2);
Xe3  = Xe(3);
omg2 = omg(2);
omg3 = omg(3);
q1   = q(1);
q2   = q(2);
q3   = q(3);
q4   = q(4);