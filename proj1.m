clc;
clear all;
close all;

% Want: thrust, specific thrust, thurst specific fuel consumption

% Table 1: Parameters for cruise (index 1) and SLS (index 2)

h = [10.67 * 10^3 0]; % m
T_0_static = [218.8 288.15]; % K
P_0_static = [0.239 1.014] * 10^5; % Pa
Ma = [0.78 0];
P_ratio_overall = [32 28];
P_ratio_fan = [1.55 1.52];
T_04 = [1450 1650]; % K
m_dot = [110 265]; % kg/s
BR = [10 10];

% Table 2: Engine Data

P_02_over_00 = [0.998 1.00];
eta_fan = [0.95 0.95];
eta_comp = [0.89 0.89];
eta_turb = [0.90 0.90];
eta_nozz = [0.95 0.95];
P_04_over_03 = [0.95 0.95];

% Part 2
k = 1.4;
c_p = 1005; % J / kg * K

T_0_stag = 
P_0_stag = 



