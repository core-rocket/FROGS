function [] = FROGSthrust
% FROGS
% ver1.8 (190807edited)
% for NSE15th
% Thrust data 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global thrust mm0
mm0 = 2.476;    %点火時下段エンジン質量[kg]
Thrustdata = readmatrix('Thrustdata0609.csv');
thrust = Thrustdata(:,2);