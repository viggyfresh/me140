clc;
clear all;
close all;

%Raw data
%Ignoring last two data points due to experiment malfunction
rpm = 70500;
Tm2 = 20.9066 + 273.15; %Cross flow
Tm3 = 171.9082 + 273.15; %Cross flow
Tm4 = 597.3990 + 273.15; %Cross flow
Tm5 = 556.9158 + 273.15; %Axial flow
Tm8 = 514.3652 + 273.15; %Cross flow
dp2 = 1.6022 * 10^3; %Differential
pt3 = 132.8239 * 10^3; %Stagnation
p4 =  128.2483 * 10^3; %Static
pt5 = 15.7876 * 10^3; %Stagnation
pt8 = 11.2959 * 10^3; %Stagnation
m_dot_fuel = 0.0032; %kg/s
thrust = [3.7 3.9 5.1 6 6.5 7.9] * 4.4482216; %N

Po2 = 101.3 * 10^3; % Pa

%Convert pressures from gauge to absolute
pt3 = pt3 + Po2;
p4 = p4 + Po2;
pt5 = pt5 + Po2;
pt8 = pt8 + Po2;

%Given/known information
A1 = 27.3 * 0.00064516;
A2 = 6.4 * 0.00064516;
A3 = 9 * 0.00064516;
A4 = 7.2 * 0.00064516;
A5 = 4.7 * 0.00064516;
A8 = 3.87 * 0.00064516;
A8 = A8 * (0.5:0.25:1.5); 
RF_c = 0.68;
RF_a = 0.86;

%Time to actually find air m_dot, Ma, U, and rho at state 2
%Assumption - since Ma will be small, T2 = T2_measured ~= T2_actual
[~, ~, k, R] = sp_heats(Tm2, 'air');
Po2_over_P = Po2 ./ (Po2 - dp2);
Ma_2 = sqrt((Po2_over_P.^((k - 1) ./ k) - 1) .* (2 ./ (k - 1)));
U_2 = sqrt(k .* R .* Tm2) .* Ma_2;
rho_2 = (Po2 - dp2) ./ (R .* Tm2);
m_dot = rho_2 .* U_2 .* A2;

%Calculate air-fuel ratio
af = m_dot ./ m_dot_fuel;

%Define LHV
LHV = (42800 * 10^3) * 170.145/1000; %converted to J/mol

MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2O = 18.016;
MM.CO2 = 44.01;
MM.JetA = 170.145;

%stochiometric air fuel and equivalence ratio
AF_s = (17.85 * MM.O2 + 17.85*(79/21) * MM.N2) / (12.3 * MM.C + 22.2 * MM.H);
phi = AF_s ./ af;

for i=1:length(A8)
    [Ma2(i), To2(i), T2(i), Po2_ratio(i)] = ...
        zachStuart(Tm2, Po2, m_dot, A2, RF_c, 'air');
    [Ma3(i), To3(i), T3(i), Po3_ratio(i)] = ...
        zachStuart(Tm3, pt3, m_dot, A3, RF_c, 'air');
    %assume static = stagnation pressure at station 4 due to low Ma
    [Ma4(i), To4(i), T4(i), Po4_ratio(i)] = ...
        viggyFresh(Tm4, p4, m_dot + m_dot_fuel, A4, RF_c, phi, MM);
    [Ma5(i), To5(i), T5(i), Po5_ratio(i)] = ...
        viggyFresh(Tm5, pt5, m_dot + m_dot_fuel, A5, RF_a, phi, MM);
    [Ma8(i), To8(i), T8(i), Po8_ratio(i)] = ...
        zachStuart(Tm8, pt8, m_dot + m_dot_fuel, A8(i), RF_c, 'air');
end

%Find station 1 values
Po1 = Po2; 
To1 = To2;
for i = 1:length(A8)
    [Ma1(i), T1(i), Po1_ratio(i)] = richieTran(To1(i), Po1, m_dot, A1);
end

%Calculate speed of sound at each station
[~, ~, gamma1, ~] = sp_heats(T1, 'air');
[~, ~, gamma2, ~] = sp_heats(T2, 'air');
[~, ~, gamma3, ~] = sp_heats(T3, 'air');
[~, ~, gamma4, ~] = sp_heats_JetA(T4, phi, MM);
[~, ~, gamma5, ~] = sp_heats_JetA(T5, phi, MM);
[~, ~, gamma8, ~] = sp_heats(T8, 'air');

U1 = Ma1 .* sqrt(gamma1 .* R .* T1);
U2 = Ma2 .* sqrt(gamma2 .* R .* T2);
U3 = Ma3 .* sqrt(gamma3 .* R .* T3);
U4 = Ma4 .* sqrt(gamma4 .* R .* T4);
U5 = Ma5 .* sqrt(gamma5 .* R .* T5);
U8 = Ma8 .* sqrt(gamma8 .* R .* T8);

%Compute all static pressures from stagnation ratios
stag_two = ones(1, length(A8)) * Po2;
P1 = stag_two ./ Po1_ratio;
P2 = stag_two ./ Po2_ratio;
P3 = pt3 ./ Po3_ratio;
P4 = p4./ Po4_ratio;
P5 = pt5 ./ Po5_ratio;
P8 = pt8 ./ Po8_ratio;

%Calculate thrust terms - CV from state 0 to state 8
Ft_calc = ((m_dot + m_dot_fuel) .* U8);
thrust_sp = Ft_calc ./ m_dot;
TSFC = m_dot_fuel ./ Ft_calc;

%Marker size var
markerSize = 10;

%Convert rpm tp krmp
krpm = rpm ./ 1000;

%Find Q_dot into system and work out of turbine
lhv = 42800 * 10^3; %J/kg
Q_dot = m_dot_fuel .* lhv;
W_net = (m_dot + m_dot_fuel) .* (U8 .^ 2)  ./ 2;
eta_therm = W_net ./ Q_dot;

figure;
plot(A8, Ft_calc, 'marker','o','MarkerSize',markerSize);
xlabel('Exit Nozzle Area (m^2)');
ylabel('Calculated Thrust (N)');
title('Calculated Thrust vs. Exit Nozzle Area');
set(gcf,'color','w');
set(gca, 'XTickLabel', num2str(get(gca,'XTick')', '%f'));

figure;
plot(A8, thrust_sp, 'marker','o','MarkerSize',markerSize);
xlabel('Exit Nozzle Area (m^2)');
ylabel('Specific Thrust (N*s/kg)');
title('Specific Thrust vs. Exit Nozzle Area');
set(gcf,'color','w');
set(gca, 'XTickLabel', num2str(get(gca,'XTick')', '%f'));

plotfixer;