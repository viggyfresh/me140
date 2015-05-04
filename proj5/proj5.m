clear all;
clc;
close all;

%% Part A1

% measured values (that is, TA data measured values!)
load = [1 2 3 4 5];
resistors = [0 1 2 3 4];
V_stack = [14.82 13.21 12.15 10.46 9.21]; %volts
I_stack = [5.55 13.43 19.84 30.3 38.84]; %amps
V_load = [14.77 13.15 12.04 10.38 8.98]; %volts
I_load = [0 6.5 11.74 19.95 25.3]; %amps

H2_flow = [2.5 4.3 6 10 14] .* (1/3600) * (.3048^3) * (.0899); %SCFH to kg/s
air_flow = [.6 .75 .75 1.1 1.3] .* (1/60) * (.3048^3) * (1.293); %SCFM to kg/s

T_air_in_stack = [49.2 48.7 48.7 49.9 52.3]; %celsius
T_air_out_stack = [46.6 46.4 46.4 47.7 50];
T_water_reservoir = [51 50.5 50.2 50.7 51.4];
T_water_in_stack = [51.1 50.4 50.3 50.9 51.3];
T_water_before_HeatExchange = [50.5 50.2 50.5 51.8 53.5];
T_stack = [48.2 47 47.1 47.6 49.4];
P_air_in = [.75 1 1.2 1.6 2]; %these are GAUGE!!!
P_H2_in = [1 1.1 1.1 1.1 1.1]; %GAUGE


%Plots for Part 1
P_load = I_load .* V_load; %watts
P_stack = I_stack .* V_stack; %watts
P_accessory = P_stack - P_load; %watts

figure
plot(P_load, I_load, P_load, I_stack);
xlabel('Power to Resistor Bank (W)');
ylabel('Current (A)');
title('Load and Stack Currents vs. Load Power');
legend('Load current', 'Stack current','Location','northwest');
set(gcf, 'color', 'w');
plotfixer; %% this must be added 

figure 
plot(P_load, V_load, P_load, V_stack);
xlabel('Power to Resistor Bank (W)');
ylabel('Potential/Voltage (V)');
title('Load and Stack Potentials vs. Load Power');
legend('Load potential', 'Stack potential','Location','southwest');
set(gcf, 'color', 'w');
plotfixer; 

figure
plot(P_load, P_stack, P_load, P_accessory); % need to put where net power is zero
xlabel('Power to Resistor Bank (W)');
ylabel('Power (W)');
title('Stack and Accessory Power vs. Load Power');
legend('Stack Power', 'Accessory Power');
set(gcf, 'color', 'w');
plotfixer;

figure 
plot(P_load, H2_flow, P_load, air_flow);
xlabel('Power to Resistor Bank (W)');
ylabel('Mass flow rate (kg/s)');
title('Mass flow rate of Hydrogen and Air vs. Load Power');
legend('Hydrogen gas', 'Air','Location','northwest');
set(gcf, 'color', 'w');
plotfixer; 

%% Part A2



% Lambda vs. I_load
figure 
plot(I_load, lambda)
xlabel('Excess-air coefficient');
ylabel('Load Current [volts]')
plotfixer 

% n_1 vs. I_load


% n_2 vs. I_load 

