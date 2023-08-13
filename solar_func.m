function [Iout,FF] = solar_func(Irr,Temp)
close all;
%% parameter
e = 1.6e-19;
kb = 1.38e-23;
T = 300; % temperature in kelvin
I0 = 25e-9; %reverse saturation current(A)
K = 2e-5; % for Si solar cell
Iph = K*Irr; %Photocurrent of the solar cell(A)
Eg = 1.14; % Bandgap of Si
n = 1; %Ideality factor
Rs = 10; %series resistance of the equivalent model
Rp = 1e6; %parallel resistance of the equivalent model
%% Temperature Effect
T_new = Temp;
I0_old = I0;
I0=((T_new^3)*exp(Eg./(kb*T_new/e)))*I0_old/((T^3)*exp(
Eg/(kb*T/e))); %reverse saturation current(A)
%% Calculation of current, power(considering Rp and Rs)
V=0:0.0001:0.4;
I_total = zeros(1,length(V));
for i = 1:length(V)
fcn = @(I) -I - Iph + I0*(exp(e*(V(i)-
I*Rs)/(n*kb*T))-1) + (V(i)-I*Rs)/Rp;
I = fzero(fcn,Iph);
I_total(i)= I;
end
Power = (-I_total.*V);
%% Load Line
R = 13.9476;
err = (-V/R-I_total);32
index = find(abs(err) == min(abs(err)));
%% I-V Curve Plot
figure
plot(V,I_total*1e3,'Linewidth',2)
xlabel('Voltage, V(V)', 'FontWeight','bold')
ylabel('Current,I_{total}(mA)', 'FontWeight','bold')
grid on;
hold on
%line([V(1), V(end)], [0, 0], 'Color',
[0,0,0],'LineStyle','-.','linewidth',2);
plot(V,(-V/R)*1e3);
plot(V(index),I_total(index)*1e3,'ro');
title('I-V Characteristics of Solar Cell with Load
Line')
legend({'Photodiode I-V', 'Load Line', 'Operating
Point'}, 'FontWeight','bold')
Iout = I_total(index);
% Vout = V(index);
% Pout = (-Iout)*Vout;
% figure
% plot(V,Power)
% xlabel('Voltage, V(V)')
% ylabel('Power(W)')
%% Fill Factor calculation
Isc = Iph;
Pmax = max(Power);
index = find(min(abs(I_total)) == abs(I_total));
Voc = V(index);
FF = Pmax/(Isc*Voc);