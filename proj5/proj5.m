clear all
clc
close all

%measured values (that is, TA data measured values!)
load = [1 2 3 4 5];
resistors = [0 1 2 3 4];
V_stack = [14.82 13.21 12.15 10.46 9.21]; %volts
I_stack = [5.55 13.43 19.84 30.3 38.84]; %amps
V_load = [14.77 13.15 12.04 10.38 8.98]; %volts
I_load = [0 6.5 11.74 19.95 25.3]; %amps
H2_flow = [29 29 29 29 29]; %SCFH
air_flow = [.6 .75 .75 1.1 1.3]; %SCFH
T_air_in_stack = [49.2 48.7 48.7 49.9 52.3]; %celsius
T_air_out_stack = [46.6 46.4 46.4 47.7 50];
T_water_reservoir = [51 50.5 50.2 50.7 51.4];
T_water_in_stack = [51.1 50.4 50.3 50.9 51.3];
T_water_before_HeatExchange = [50.5 50.2 50.5 51.8 53.5];
T_stack = [48.2 47 47.1 47.6 49.4];
P_air_in = [75 1 1.2 1.6 2];
P_H2_in = [1 1.1 1.1 1.1 1.1];

