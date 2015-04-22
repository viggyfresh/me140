%Part 4, Varying paramters
clc;
clear all;
close all;

% varying phi
phi = [0.0911	0.1821	0.2732	0.3642	0.7285	0.9106];
Thrust_phi = [50.4632	51.1187	51.78	52.4464	55.1681	56.561];
markerSize = 10;
plot(phi, Thrust_phi, 'marker', 'o', 'MarkerSize', markerSize);
xlabel('Equivalence Ratio (phi))');
ylabel('Thrust (N)');
title('Thrust vs. Equivalence Ratio');
set(gcf,'color','w');

%Raw data at max spool speed
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
thrust = 7.9 * 4.4482216; %N

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

% Adiabatic flame temp is about 900K
To4 = [700 750 800 850 900];
eta_turb = 0.2879;
eta_nozz = 0.5848;

for i=1:length(To4)
    To4(i)
    [Ma2(i), To2(i), T2(i), Po2_ratio(i)] = ...
        zachStuart(Tm2, Po2, m_dot, A2, RF_c, 'air');
    [Ma3(i), To3(i), T3(i), Po3_ratio(i)] = ...
        zachStuart(Tm3, pt3, m_dot, A3, RF_c, 'air');
    
    [Ma4(i), T4(i), Po5_ratio(i)] = richieTran(To4(i), p4, m_dot + m_dot_fuel, A4);
   
    To5s(i) = turb_Ts(To4(i), pt5 / p4, 1, 'JetA', phi, MM);
    To5(i) = var_cp_turb(To4(i), To5s(i), eta_turb);
    [Ma5(i), T5(i), Po5_ratio(i)] = richieTran(To5(i), pt5, m_dot + m_dot_fuel, A5);
    
    To8s(i) = var_cp_neg(To5(i), pt8 / pt5);
    To8(i) = var_cp_nozz(To5(i), To8s(i), eta_nozz);
    [Ma8(i), T8(i), Po8_ratio(i)] = richieTran(To8(i), pt8, m_dot + m_dot_fuel, A8);

    
end

%Calculate speed of sound at each station
[~, ~, gamma2, ~] = sp_heats(T2, 'air');
[~, ~, gamma3, ~] = sp_heats(T3, 'air');
[~, ~, gamma4, ~] = sp_heats_JetA(T4, phi, MM);
[~, ~, gamma5, ~] = sp_heats_JetA(T5, phi, MM);
[~, ~, gamma8, ~] = sp_heats_JetA(T8, phi, MM);


U2 = Ma2 .* sqrt(gamma2 .* R .* T2);
U3 = Ma3 .* sqrt(gamma3 .* R .* T3);
U4 = Ma4 .* sqrt(gamma4 .* R .* T4);
U5 = Ma5 .* sqrt(gamma5 .* R .* T5);
U8 = Ma8 .* sqrt(gamma8 .* R .* T8);

%Compute all static pressures from stagnation ratios
stag_two = ones(1, length(rpm)) * Po2;
P2 = stag_two ./ Po2_ratio;
P3 = pt3 ./ Po3_ratio;
%P4 = p4./ Po4_ratio;
P5 = pt5 ./ Po5_ratio;
P8 = pt8 ./ Po8_ratio;

%Calculate thrust terms - CV from state 0 to state 8
Ft_calc = ((m_dot + m_dot_fuel) .* U8)
thrust_sp = Ft_calc ./ m_dot
TSFC = m_dot_fuel ./ Ft_calc

%Marker size var
markerSize = 10;

%Convert rpm tp krmp
krpm = rpm ./ 1000;

figure;
plot(To4, Ft_calc, 'marker', 'o', 'MarkerSize', markerSize);
xlabel('Turbine Inlet Temperature To4 (K)');
ylabel('Thrust (N)');
title('Thrust vs. Turbine Inlet Temperature');
set(gcf,'color','w');

%Find Q_dot into system and work out of turbine
lhv = 42800 * 10^3; %J/kg
Q_dot = m_dot_fuel .* lhv;
W_net = (m_dot + m_dot_fuel) .* (U8 .^ 2)  ./ 2;
eta_therm = W_net ./ Q_dot

figure;
plot(To4, eta_therm * 100, 'marker', 'o', 'MarkerSize', markerSize);
xlabel('Thermal Efficiency vs. Inlet Temperature To4 (K)');
ylabel('Thermal Efficiency (%)');
title('Thermal Efficiency vs. Turbine Inlet Temperature');
set(gcf,'color','w');

plotfixer;