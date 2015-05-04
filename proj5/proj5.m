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
P_air_in = TA_data(:,9)' + P_atm;  % Converted to absolute
P_H2_in = TA_data(:,11)' + P_atm;  % Converted to absolute

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

% lambda vs. I_load
MM.O2 = 32;
MM.N2 = 28.02;
MM.H = 1.008;
MM.H2 = 2 * MM.H;
MM.air = 28.97;

H2_flow_m_s = H2_flow ./ MM.H2 .* 1000;  %mol/sec
Air_flow_m_s = air_flow ./MM.air .* 1000; %mol/sec
% divide by 4.76 to account for number of moles of air
lambda = (2*Air_flow_m_s./H2_flow_m_s) ./ 4.76; % check with TA

figure;
plot(I_load, lambda);
xlabel('Load Current [volts]');
ylabel('\lambda');
title('\lambda vs. Load Current');
set(gcf, 'color', 'w');
plotfixer;

% n_1 vs. I_load
for i=1:length(T_stack)
    alpha(i) = relHumidity(T_stack(i),lambda(i));
    eta_1_LHV(i) = lucio(T_stack(i),P_air_in(i),P_H2_in(i),alpha(i),lambda(i))
end


figure;
plot(I_load, eta_1_LHV);
xlabel('Load Current [volts]');
ylabel('\eta (First Law)');
title('First Law Efficiency vs. Load Current');
set(gcf, 'color', 'w');
plotfixer;


% n_2 vs. I_load 

