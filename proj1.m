clc;
clear all;
close all;

% Want: thrust, specific thrust, thurst specific fuel consumption

% Table 1: Parameters for cruise (index 1) and SLS (index 2)

vars.h = [10.67 * 10^3 0]; % m
vars.T_0_static = [218.8 288.15]; % K
vars.P_0_static = [0.239 1.014] * 10^5; % Pa
vars.Ma = [0.78 0];
vars.P_ratio_overall = [32 28];
vars.P_ratio_fan = [1.55 1.52];
vars.T_04 = [1450 1650]; % K
vars.m_dot = [110 265]; % kg/s
vars.BR = [10 10];

% Table 2: Engine Data

vars.P_02_over_00 = [0.998 1.00];
vars.eta_fan = [0.95 0.95];
vars.eta_comp = [0.89 0.89];
vars.eta_turb = [0.90 0.90];
vars.eta_nozz = [0.95 0.95];
vars.P_04_over_03 = [0.95 0.95];

% Part 2
vars.k = 1.4;
vars.c_p = 1005; % J / kg * K






