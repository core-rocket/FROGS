function [] = FROGSthrust
% FROGS
% ver1.8 (190807edited)
% for NSE15th
% Thrust data 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global thrust mm0
mm0 = 2.476;    %�_�Ύ����i�G���W������[kg]
Thrustdata = readmatrix('Thrustdata0609.csv');
thrust = Thrustdata(:,2);