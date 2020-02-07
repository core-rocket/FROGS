function [] = FROGSthrust
% FROGS
% ver1.8 (190807edited)
% for NSE15th
% Thrust data 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global thrust mm0
mm0 = 6.787;    %点火時下段エンジン質量[kg]
Thrustdata = readmatrix('thrust/Thrustdata_sim_LIATRIS1100.csv');
thrust = Thrustdata(:,2);
