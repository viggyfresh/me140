%Part 1, find a new Cp

clc;
clear all;

%Raw data
%Ignoring last two data points due to experiment malfunction
rpm = [46500 49300 55000 60000 65100 70500];
Tm2 = [21.559 20.9850 21.1541 21.2047 20.8731 20.9066] + 273.15; %Cross flow
Tm3 = [68.5792 74.8547 86.1913 117.9889 129.9338 148.9507 171.9082]...
       + 273.15; %Cross flow
Tm4 = [503.7978 501.4857 509.7827 535.9689 547.0686 597.3990]...
       + 273.15; %Cross flow
Tm5 = [484.6335 474.2929 482.9551 507.9552 513.6135 556.9158]...
       + 273.15; %Axial flow
Tm8 = [488.5296 484.472 486.4055 494.6287 500.1043	514.3652]...
       + 273.15; %Cross flow
Tm_oil = [45.6508 49.6376 57.1879 64.9292 70.94 75.7677] + 273.15;
dp2 = [0.5141 0.6161 0.8152 1.0346 1.3044 1.6022] * 10^3; %Differential
pt3 = [50.4636 57.6214 73.4094 89.7863 109.6062 132.8239] * 10^3; %Stagnation
p4 = [46.6808 53.9519 69.468 85.8223 105.2571 128.2483] * 10^3; %Static
pt5 = [5.4178 6.1699 7.9933 10.2856 12.6753 15.7876] * 10^3; %Stagnation
pt8 = [3.5038 4.2604 5.9569 7.7169 9.3909 11.2959] * 10^3; %Stagnation
m_dot_fuel = [0.0021 0.0023 0.0025 0.0027 0.0029 0.0032]; %kg/s
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
RF_c = 0.68;
RF_a = 0.86;

R = 286.9;

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

for i=1:length(rpm)
    [Ma2(i), To2(i), T2(i), Po2_ratio(i)] = ...
        zachStuart(Tm2(i), Po2, m_dot(i), A2, RF_c);
    [Ma3(i), To3(i), T3(i), Po3_ratio(i)] = ...
        zachStuart(Tm3(i), pt3(i), m_dot(i), A3, RF_c);
end

%Define LHV
LHV = (42800 * 10^3) * 170.145/1000; %converted to J/mol

%Chemistry, NOT FOR PART 1

MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2O = 18.016;
MM.CO2 = 44.01;
MM.LHV = 170;

%stochiometric air fuel and equivalence ratio
AF_s = (17.85 * MM.O2 + 17.85*(79/21) * MM.N2) / (12.3 * MM.C + 22.2 * MM.H);
phi = AF_s ./ af;

%%%%% Find Temperature across combustor: To4 %%%%%%%%
hf.H2O = -241820; %for vapor, in J/mol 
hf.CO2 = -393520; %in J/mol
hf.JetA = 11.1 * hf.H2O + 12.3 * hf.CO2 + LHV;


%To4 = tempCalc_combustor(hf);
for i=1:length(rpm)
To4(i) = combustor(MM, phi(i), To3(i));
end
To4

% PART 2 %
[hf_mol, hf_kg] = heatOfFormation();
T_a = flameTemp(0.319, 'JetA', hf_mol, MM)
