clc;
clear all;
load('bradyfire1.mat');

rho = 950 * 3.6127292 * 10^-5;  %lb / in^3;


% Regression Stuff

%% TA Fire
load('bradycheated.mat');
i1 = start_index;
i2 = final_index;
t_2 = time(i2) - time(i1); % secs
mfuel_2 = mfuel / 10^3 * 2.205; % lbm;
m_dot_O2_2 = mean(m_dot_O2(i1:i2)) * 2.205;  % lbm/s

% dimensions of nozzle
D_center = 0.605; % inches
L_2 = 5.375; % inches
A_port_2 = pi * D_center^2 / 4; 

% calculate r_dot and Go for fitting
delta_V_2 = mfuel_2 / rho %change in Volume
delta_r_2 = sqrt(delta_V_2 / (pi * L_2))
r_dot_2 = delta_r_2 / t_2   % in/sec

Go_2 = m_dot_O2_2 / A_port_2 % lbm/(s*in^2)

%% data point 1: First Fire
load('bradyfire1.mat');
i1 = start_index;
i2 = final_index;
t_1 = time(i2) - time(i1); % secs
mfuel_1 = mfuel / 10^3 * 2.205; % lbm
m_dot_O2_1 = mean(m_dot_O2(i1:i2)) * 2.205; %lbm/s

% dimensions of nozzle
D_center = 0.5; % inches
D_outer = 0.30; % inches
L_1 = 5.375; % inches
A_port_1 = pi * D_center^2 / 4 + 8 * pi * D_outer^2 / 4; % cross_sectional area

% calculate r_dot and Go for fitting
delta_V_1 = mfuel_1 / rho; %change in Volume
delta_r_1 = sqrt(delta_V_1 / (pi * L_1))
r_dot_1 = delta_r_1 / t_1 % in/s

Go_1 = m_dot_O2_1 / A_port_1 % lbm/(s*in^2)



%% Fit regression line
% x = [Go_1 Go_2];
%Note: This is using data from "I LOVE DARREN" right now
x = [Go_1 0.1844 0.2385];
y = [r_dot_1 0.1380 0.2076];
%y = [r_dot_1 r_dot_2];

p = polyfit(log(x), log(y), 1);
n = p(1)
a = exp(p(2))


