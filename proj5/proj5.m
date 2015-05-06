clear all;
clc;
close all;

%% Part A1

% measured values (that is, TA data measured values!)
P_atm = 101325; % Pa
TA_data = xlsread('TA_data.xlsx');
I_load = TA_data(:,1)'; %amps
I_stack = TA_data(:,2)'; 
V_load = TA_data(:,3)'; %volts
V_stack = TA_data(:,4)';
H2_flow = TA_data(:,5)' .* (1/3600) .* (.3048^3) .* (.0899); %SCFH to kg/s
air_flow = TA_data(:,6)' .* (1/60) .* (.3048^3) .* (1.293); %SCFM to kg/s
T_air_in_stack = TA_data(:,12)' + 273.15; % K
T_air_out_stack = TA_data(:,13)' + 273.15; % K
T_water_reservoir = TA_data(:,14)' + 273.15; % K
T_water_in_stack = TA_data(:,15)' + 273.15; % K
T_water_before_HeatExchange = TA_data(:,16)' + 273.15; % K
T_stack = TA_data(:,17)' + 273.15; % K 
P_air_in = TA_data(:,9)' + P_atm;  % Converted to absolute (Pa)
P_H2_in = TA_data(:,11)' + P_atm;  % Converted to absolute (Pa)

%Plots for Part 1
P_load = I_load .* V_load; %watts
P_stack = I_stack .* V_stack; %watts
P_accessory = P_stack - P_load; %watts

figure;
plot(P_load, I_load, P_load, I_stack);
xlabel('Power to Resistor Bank (W)');
ylabel('Current (A)');
title('Load and Stack Currents vs. Load Power');
legend('Load current', 'Stack current','Location','northwest');
set(gcf, 'color', 'w');
plotfixer; %% this must be added 

figure;
plot(P_load, V_load, P_load, V_stack);
xlabel('Power to Resistor Bank (W)');
ylabel('Potential/Voltage (V)');
title('Load and Stack Potentials vs. Load Power');
legend('Load potential', 'Stack potential','Location','southwest');
set(gcf, 'color', 'w');
plotfixer; 

figure;
plot(P_load, P_stack, P_load, P_accessory); % need to put where net power is zero
xlabel('Power to Resistor Bank (W)');
ylabel('Power (W)');
title('Stack and Accessory Power vs. Load Power');
legend('Stack Power', 'Accessory Power');
set(gcf, 'color', 'w');
plotfixer;

figure;
plot(P_load, H2_flow, P_load, air_flow);
xlabel('Power to Resistor Bank (W)');
ylabel('Mass flow rate (kg/s)');
title('Mass flow rate of Hydrogen and Air vs. Load Power');
legend('Hydrogen gas', 'Air','Location','northwest');
set(gcf, 'color', 'w');
plotfixer; 

%% Part A2

% lambda vs. P_load
MM.O2 = 32;
MM.N2 = 28.02;
MM.H = 1.008;
MM.H2 = 2 * MM.H;
MM.air = 28.97;

H2_flow_mol_s = H2_flow ./ MM.H2 .* 1000;  %mol/sec
Air_flow_mol_s = air_flow ./MM.air .* 1000; %mol/sec
% divide by 4.76 to account for number of moles of air
lambda = (2*Air_flow_mol_s./H2_flow_mol_s) ./ 4.76; % check with TA

figure;
plot(P_load, lambda);
xlabel('Power to Resistor Bank (W)');
ylabel('\lambda');
title('\lambda vs. Load Power');
set(gcf, 'color', 'w');
plotfixer;
% n_1 vs. P_load
for i=1:length(T_stack)
    alpha(i) = john(T_stack(i), lambda(i), P_air_in(i));
    [deltaG_rxn(i)] = lucio(T_stack(i),P_air_in(i),P_H2_in(i),alpha(i),lambda(i));
end


deltaG = deltaG_rxn .* H2_flow_mol_s; % per mol H2 basis
LHV = 120 * 10^6; 

eta_1_load = P_load ./ (LHV.*H2_flow);
eta_2_load = P_load ./ (-deltaG);

eta_1_stack = P_stack ./ (LHV.*H2_flow);
eta_2_stack = P_stack ./ (-deltaG);

P_lost_load = -deltaG - P_load;
P_lost_stack = -deltaG - P_stack;

figure;
plot(P_load, eta_1_load * 100, P_load, eta_1_stack * 100);
xlabel('Power (W)');
ylabel('\eta_1');
title('First Law Efficiency vs. Load');
legend('Load', 'Stack', 'Location', 'Southeast');
set(gcf, 'color', 'w');
plotfixer;

figure;
plot (P_stack, eta_2_load * 100, P_stack, eta_2_stack * 100);
xlabel('Power (W)');
ylabel('\eta_2');
title('Second Law Efficiency vs. Load');
legend('Load','Stack');
set(gcf, 'color', 'w');
plotfixer;

figure;
plot(P_load, P_lost_load, P_load, P_lost_stack);
xlabel('Power (W)');
ylabel('Lost Power (W)');
title('Lost Power vs. Load Power');
legend('Load', 'Stack');
set(gcf, 'color', 'w');
plotfixer;

%% Part B1
T_range = (25+273):10:(1200+273);
R = 8.3144621; %universal gas constant 

% only supposed to plot 10^-3 < Kp < 10^3 
% maybe split into two for loops and insert if statements?
% also could look into refining plot function
for i = 1:length(T_range)
    deltaG_reform(i) = lucio_smr(T_range(i));
    deltaG_wgs(i) = lucio_wgs(T_range(i));
    Kp_reform(i) = exp(-deltaG_reform(i) ./ (R*T_range(i)));
    Kp_wgs(i) = exp(-deltaG_wgs(i) ./ (R*T_range(i)));
    i
end 

figure;
semilogy(T_range, Kp_reform, T_range, Kp_wgs);
xlabel('Temperature [K]');
ylabel('Equilibrium Constant (Kp)');
title('Temp vs. Kp');
set(gcf, 'color', 'w');
legend('K_p SMR','K_p WGS','Location','Northwest');
plotfixer;


