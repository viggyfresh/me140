clc
close all
clear all

%Raw data
rpm = [46500 49300 55000 60000 65100 70500]; %Ignoring last two data points due to experiment malfunction
Tm2 = [21.559 20.9850 21.1541 21.2047 20.8731 20.9066] + 273.15; %Cross flow
Tm3 = [68.5792 74.8547 86.1913 117.9889 129.9338 148.9507 171.9082] + 273.15; %Cross flow
Tm4 = [503.7978 501.4857 509.7827 535.9689 547.0686 597.3990] + 273.15; %Cross flow
Tm5 = [484.6335 474.2929 482.9551 507.9552 513.6135 556.9158] + 273.15; %Axial flow
Tm8 = [488.5296 484.472 486.4055 494.6287 500.1043	514.3652] + 273.15; %Cross flow
Tm_oil = [45.6508 49.6376 57.1879 64.9292 70.94 75.7677] + 273.15;
dp2 = [0.5141 0.6161 0.8152 1.0346 1.3044 1.6022] * 10^3; %Stagnation - static (differential)
pt3 = [50.4636 57.6214 73.4094 89.7863 109.6062 132.8239] * 10^3; %Stagnation
p4 = [46.6808 53.9519 69.468 85.8223 105.2571 128.2483] * 10^3; %Static
pt5 = [5.4178 6.1699 7.9933 10.2856 12.6753 15.7876] * 10^3; %Stagnation
pt8 = [3.5038 4.2604 5.9569 7.7169 9.3909 11.2959] * 10^3; %Stagnation
fuelFlow = [0.0021 0.0023 0.0025 0.0027 0.0029 0.0032]; %kg/s
thrust = [3.7 3.9 5.1 6 6.5 7.9] * 4.4482216; %N

Po2 = 101.3 * 10^3; % Pa

%Given/known information
A1 = 27.3 * 0.00064516;
A2 = 6.4 * 0.00064516;
A3 = 9 * 0.00064516;
A4 = 7.2 * 0.00064516;
A5 = 4.7 * 0.00064516;
A8 = 3.87 * 0.00064516;
RF_c = 0.68;
RF_a = 0.86;

R = 286.9;

%Does the chemistries
M_O2 = 32; %g/mol
M_N2 = 28.013;
M_C = 12.011;
M_H = 1.008;

%Calculate air/fuel ratio
AF_ratio = 18.5 * ((M_O2 + (79/21) * M_N2)) / (12 * M_C + 26 * M_H);

%Find air mass flow - not sure why
airFlow = fuelFlow*AF_ratio;

%Time to actually find air m_dot, Ma, U, and rho at state 2
%Assumption - since Ma will be small, T2 = T2_measured ~= T2_actual

[~, ~, k, R] = sp_heats(T2);
Po2_over_P = Po2 ./ (Po2 - dp2);
Ma_2 = sqrt((Po2_over_P.^((k - 1) ./ k) - 1) .* (2 ./ (k - 1)));
U_2 = sqrt(k .* R .* T2) .* Ma_2;
rho_2 = (Po2 - dp2) ./ (R .* T2);
m_dot = rho_2 .* U_2 .* A2;

%Find mach number in order to find static and stagnation temperature values
for i=1:6
    Ma2(i), To2(i), T2(i) = stationAnalysis(Tm2(i), Po2, m_dot(i), A2, RF_c.68)
    Ma3(i), To3(i), T3(i) = stationAnalysis(Tm3(i), pt3(i), m_dot(i), A3, RF_c)
    %assume static = stagnation pressure at station 4 due to low Ma
    Ma4(i), To4(i), T4(i) = stationAnalysis(Tm4(i), p4(i), m_dot(i), A4, RF_c)
    Ma5(i), To5(i), T5(i) = stationAnalysis(Tm5(i), pt5(i), m_dot(i), A5, RF_a)
    Ma8(i), To8(i), T8(i) = stationAnalysis(Tm8(i), pt8(i), m_dot(i), A8, RF_c)
end
