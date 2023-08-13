clc;
close all;
clear all;
%% Parameter
Temp = 290; %temperature in kelvin
Irr = 500; %Irradiance(Wm-2)
%% Solar cell
[Iout_sc,FF] = solar_func(Irr,Temp);
fprintf("Output current from solar cell=%.2f mA, Fill
factor=%.2f\n", Iout_sc*1e3, FF);
%% Laser
I = -Iout_sc;
[lambda_in,Intensity,Rd] = laser_func(I);
fprintf("Output intensity from laser=%.2f MW per square
meter at wavvelength %.2f micrometer\n", Intensity/1e6,
lambda_in*1e6);
%% Photodetector
%[Iout_pd] =
photodetector_func(lambda_in,Pout_laser,Temp)
[Iout,SNR_db] =
pin_photodetector_func(lambda_in,Intensity,Temp);
fprintf("Output current from detector=%.2f nA, SNR in
dB=%.2f", Iout*1e9, SNR_db);