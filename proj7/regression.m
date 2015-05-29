clc;
clear all;
load('bradyfire1.mat');

rho = 950; %kg / m^3;


% Regression Stuff

%% data point 1: First Fire
load('bradyfire1.mat');
i1 = start_index;
i2 = final_index;
t_1 = time(i2) - time(i1); %secs
mfuel_1 = mfuel / 10^3; %kg
m_dot_fuel_1 = mfuel_1 / t_1; %kg
m_dot_O2_1 = mean(m_dot_O2(i1:i2)); 

%dimensions of nozzle
D_center = 0.5 / 39.370; %meters
D_outer = 0.30  /39.370; %meters
L_1 = 5.375 / 39.370; %meters
A_port_1 = pi * D_center^2 / 4 + 8 * pi * D_outer^2 / 4; %cross_sectional area

Go_1 = m_dot_O2_1 / A_port_1;
delta_V_1 = mfuel_1 / rho; %change in Volume
delta_r_1 = sqrt(delta_V_1 / (pi * L_1));
r_dot_1 = delta_r_1 / t_1

%% TA Fire
load('bradycheated.mat');
il = start_index;
i2 = final_index;
t_2 = time(i2) - time(i1); %secs
mfuel_2 = mfuel / 10^3; %kg;
m_dot_fuel_2 = mfuel_2 / t_2; %kg/s
m_dot_O2_2 = mean(m_dot_O2(i1:i2));  %kg/s

%dimensions of nozzle
D_center = 0.605 / 39.24; %meters
L_2 = 5.375 / 39.370;
A_port_2 = pi * D_center^2 / 4;

Go_2 = m_dot_O2_2 / A_port_2;
delta_V_2 = mfuel_2 / rho %change in Volume
delta_r_2 = sqrt(delta_V_2 / (pi * L_2));
r_dot_2 = delta_r_2 / t_2

%% Fit regression line


