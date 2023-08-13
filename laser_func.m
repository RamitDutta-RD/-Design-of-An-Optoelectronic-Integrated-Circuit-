function [lambda_in,Intensity_out,Rd] = laser_func(I)
close all
%% parameter definition
N0 = 6.5e23;
a0 = 3.13e-20;
a_bar = 1.2;
lambda_0 = 1.575e-6;
b0 = 3.17e-26*1e-6;
b1 = 0;
lambda_z0 = 1.625e-6;
z0 = -2.5e-27*1e-6;
%% plotting for discrete concentration values
N = [3 2.5 2 1]*1e24;
lambda = [1.45:0.01:1.65]*1e-6;
hw = 1243./(lambda*1e9);
gp = zeros(1,length(N));
lambda_p = zeros(1,length(N));
lambda_z = zeros(1,length(N));
cn = zeros(1,length(N));
dn = zeros(1,length(N));
gm = zeros(length(lambda), length(N));
for j=1:length(lambda)
for i=1:length(N)
gp(i) = a0*(N(i)-N0) + a0*a_bar*N0*exp(-
N(i)/N0);
lambda_p(i) = lambda_0 - (b0*(N(i)-N0) +
b1*((N(i)-N0)^2));
lambda_z(i) = lambda_z0 - z0*(N(i)-N0);
cn(i) = 3*gp(i)/((lambda_z(i)-lambda_p(i))^2);34
dn(i) = 2*gp(i)/((lambda_z(i)-lambda_p(i))^3);
if lambda(j)<lambda_z(i)
gm(j,i) = cn(i)*((lambda(j) -
lambda_z(i)).^2) + ...
dn(i)*((lambda(j) - lambda_z(i)).^3);
end
end
end
figure()
for i=1:length(N)
plot(hw, gm(:,i), 'LineWidth', 2);
hold on
end
xlabel('Energy (in eV)')
ylabel('Optical gain')
legend('3e18', '2.5e18', '2e18', '1e18')
title('Optical gain vs Energy for differnet
injections')
%% plotting for continuous concentration values
N = [0.75:0.01:3]*1e24;
lambda = [1.45:0.01:1.65]*1e-6;
hw = 1243./(lambda*1e9);
gp = zeros(1,length(N));
lambda_p = zeros(1,length(N));
lambda_z = zeros(1,length(N));
cn = zeros(1,length(N));
dn = zeros(1,length(N));
gm = zeros(length(N), length(lambda));
gpeak = zeros(1,length(N));
for i=1:length(N)
for j=1:length(lambda)
gp(i) = a0*(N(i)-N0) + a0*a_bar*N0*exp(-
N(i)/N0);35
lambda_p(i) = lambda_0 - (b0*(N(i)-N0) +
b1*((N(i)-N0)^2));
lambda_z(i) = lambda_z0 - z0*(N(i)-N0);
cn(i) = 3*gp(i)/((lambda_z(i)-lambda_p(i))^2);
dn(i) = 2*gp(i)/((lambda_z(i)-lambda_p(i))^3);
if lambda(j)<lambda_z(i)
gm(i,j) = cn(i)*((lambda(j) -
lambda_z(i)).^2) + ...
dn(i)*((lambda(j) - lambda_z(i)).^3);
end
end
gpeak(i) = max(gm(i,:));
end
figure()
plot(N*1e-6,gpeak*1e-2, 'LineWidth', 2);
xlabel('Injection density (in per cubic cm)',
'FontWeight','bold')
ylabel('Peak optical gain coefficient, g_p_e_a_k (in
per cm)', 'FontWeight','bold')
title('Dependence of peak optical gain coefficient on
injection')
ylim([0 1000])
set(gca, 'Yscale', 'log')
%% parameter definition from book example for
semiconductor laser
L = 100e-6;
W = 10e-6;
d = 0.15e-6;
gamma = 2500; %loss coefficient per meter
nr = 3.491; %In0.6 Ga0.4 As0.85 P0.15
http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/optic.html
R = ((nr-1)^2)/((nr+1)^2); %reflectance
B = 2e-16; %exercise 4.3036
e = 1.6e-19;
c = 3e8;
h = 6.626e-34;
%% calculating loss
confinement = 1;
alpha_t = gamma + (1/(2*L))*log(1/(R*R)); %total loss
gth = alpha_t/confinement;
gth_arr = gth*ones(1,length(N));
% plot(N,gpeak);
% hold on
% plot(N, gth_arr);
err = gpeak - gth_arr;
tolerance = 120;
index = find(abs(err)<tolerance);
nth = N(index); %threshold electron conc
%% plot allowed modes on the same plot as optical gain
vs lambda
figure()
plot(lambda/1e-6, gm(index,:)*1e-2, 'LineWidth', 2)
%lambda(find(max(gm(index,:))))
hold on
line([lambda(1)/1e-6 lambda(end)/1e-6], [gth*1e-2
gth*1e-2], 'Color', [0 0 0], 'LineWidth', 1.5);
xlabel('Wavelength (in um)', 'FontWeight', 'bold')
ylabel('Optical gain (in per meter)', 'FontWeight',
'bold')
title('Optical gain vs wavelength at threshold electron
concentration')
max_index = find( gm(index,:) == max(gm(index,:)) );
lambda_in = lambda(max_index);
%% radiative lifetime calculation
tau_r = 1/(B*nth);
%% Threshold current density37
Jth = nth*e*d/tau_r;
Ith = W*L*Jth;
%% Photon cavity lifetime
tau_ph = nr/(c*alpha_t);
%% Output Power
%lambda_in = 1500e-9; %from figure 4.48
Pout_slope = (h*c*c*tau_ph*(1-R)/(2*e*nr*lambda_in*L));
Pout = Pout_slope*(I-Ith);
Intensity_out = Pout/(W*d);
%% Rd calculation
h = 6.626e-34;
kb = 1.38e-23;
tau_sp = 300e-9;
Ni = 1.88e17; % in m^-3 In0.6 Ga0.4 As0.85 P0.15
T = 300;
V = 1;
Ne = Ni*exp(V*e/(2*kb*T));
Rd = (2*kb*T/e)*(tau_sp/(Ne*L*W*e*d));
end