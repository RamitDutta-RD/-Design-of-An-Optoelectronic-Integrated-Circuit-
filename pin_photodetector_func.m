function [Iout,SNR_db] = pin_photodetector_func(lambda_in,Intensity,Temp)
close all;
%% parameter
e = 1.6e-19;
kb = 1.38e-23;
n = 1; %Ideality factor
T = 300; % temperature in kelvin
I0 = 25e-9; %reverse saturation current(A)
Eg = 0.784; %In0.53Ga0.47As
h = 6.626e-34;
c = 3e8;
Tr = 1; %perfect AR coating(assume)
ni = 1;% internal quantum efficiency(assume)
alpha = 4e5; %absorption coeff(in m^-1) %Figure 5.5
diameter = 0.12e-3; % in meter
%% Temperature Effect
T_new = Temp;
I0_old = I0;
I0=((T_new^3)*exp(Eg./(kb*T_new/e)))*I0_old/((T^3)*exp(
Eg/(kb*T/e))); %reverse saturation current(A)
%% Pout calculation
lambda_in = lambda_in*1e9; % from laser (in nm)
dia = 0.4e-6; % in meter
Area = (pi/4)*dia^2;
a = 10^(21/10); %considering very strong atmospheric
turbulence
dist = 5e-2;
Intensity_pd = Intensity*exp(-a*dist);
Pout = Intensity_pd*Area;
freq = c/(lambda_in*1e-9);
%% Finding Minimum W39
Width = [1.5:0.01:10]*1e-6;
Iph_W = zeros(1,length(Width));
R_W = zeros(1,length(Width));
for i=1:length(Width)
Iph(i) = e*ni*Tr*Pout*(1-exp(-
alpha*Width(i)))/(h*freq);
R_W(i) = Iph(i)/Pout;
end
plot(Width,R_W, 'LineWidth',2)
hold on
line([Width(1) Width(end)], [0.7 0.7], 'Color', [0 0
0], 'LineWidth', 1.5);
legend({'Responsivity', 'Threshold responsivity'},
'FontWeight','bold');
xlabel('Width(m)', 'FontWeight', 'bold');
ylabel('Responsivity(A/W)', 'FontWeight', 'bold');
title('Responsivity vs Width for pin Photodetector');
Rmin = 0.7; %% from the chart of the book (as we
operate near the peak)
err = R_W-Rmin; % At peak lambda Rmin = 0.7
index = find(abs(err) == min(abs(err)));
W_min = Width(index);
%% Iph calculation
W = 3e-6; % in meter
Iph = e*ni*Tr*Pout*(1-exp(-alpha*W))/(h*freq);
R = Iph/Pout;
Iph_max = e*ni*Tr*Pout/(h*freq);
%% Calculation of current, power
Vr = 1.5;
V = -2:0.0001:0;
I_total = -Iph + I0.*(exp(e*V/(n*kb*T))-1);
Power = (-I_total.*V);
%index = find(V == -Vr);
%% Load Line40
RL = 1000;
err = (-(V+Vr)/RL-I_total);
index = find(abs(err) == min(abs(err)));
%% I-V Curve Plot
figure
plot(V,I_total*1e6,'Linewidth',2)
xlabel('Voltage, V(V)')
ylabel('Current,I_{total}(uA)')
grid on;
hold on
line([V(1), V(end)], [0, 0], 'Color',
[0,0,0],'LineStyle','-.','linewidth',2);
plot(V,-((V+Vr)/RL)*1e6);
plot(V(index),I_total(index)*1e6,'ro')
title('I-V characteristics of Photodetector')
Iout = I_total(index); % in A
Vout = V(index);
%% SNR calculation
B = 1e6; %in Hz
signal_power = (Iph^2)*RL;
noise_power = [2*e*(I0+Iph)*B]*RL+4*kb*T*B;
SNR = signal_power/noise_power;
SNR_db = 10*log10(SNR);
end