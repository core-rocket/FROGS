function [] = FROGSthrust
% FROGS
% ver1.8 (190807edited)
% for NSE15th
% Thrust data 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global thrust mm0
mm0 = 3.087;    %点火時下段エンジン質量[kg]
Thrustdata = readmatrix('thrust/LIATRIS-sim.csv');
thrust = Thrustdata(:,2);
