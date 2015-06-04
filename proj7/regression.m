clc;
clear all;
load('bradyfire1.mat');

rho = 950 * 3.6127292 * 10^-5;  %lb / in^3;


% Regression Stuff

%% TA Fire
% load('bradycheated.mat');
% i1 = start_index;
% i2 = final_index;
% t_2 = time(i2) - time(i1); % secs
% mfuel_2 = mfuel / 10^3 * 2.205; % lbm;
% m_dot_O2_2 = mean(m_dot_O2(i1:i2)) * 2.205;  % lbm/s
% 
% % dimensions of nozzle
% D_center = 0.605; % inches
% L_2 = 5.375; % inches
% A_port_2 = pi * D_center^2 / 4; 
% 
% % calculate r_dot and Go for fitting
% delta_V_2 = mfuel_2 / rho %change in Volume
% delta_r_2 = sqrt(delta_V_2 / (pi * L_2))
% r_dot_2 = delta_r_2 / t_2   % in/sec
% 
% Go_2 = m_dot_O2_2 / A_port_2 % lbm/(s*in^2)

%% data point 1: First Fire
load('bradyfire1.mat');
i1 = start_index;
i2 = final_index;
t_1 = time(i2) - time(i1); % secs
mfuel_1 = mfuel / 10^3 * 2.205; % lbm
m_dot_O2_1 = mean(m_dot_O2(i1:i2)) * 2.205; %lbm/s

% dimensions of nozzle
D_center = 0.5; % inches
D_outer_1 = 0.30; % inches
L_1 = 5.375; % inches
A_port_1 = pi * D_center^2 / 4 + 8 * pi * D_outer_1^2 / 4; % cross_sectional area

% calculate r_dot and Go for fitting
delta_V_1 = mfuel_1 / rho; %change in Volume
delta_r_1 = sqrt(delta_V_1 / (pi * L_1));
r_dot_1 = delta_r_1 / t_1; % in/s

Go_1 = m_dot_O2_1 / A_port_1 % lbm/(s*in^2)

%% data point 2: 2nd Fire
load('bradyfire2.mat');
i1 = start_index;
i2 = final_index;
t_2 = time(i2) - time(i1); % secs
mfuel_2 = mfuel / 10^3 * 2.205; % lbm;
m_dot_O2_2 = mean(m_dot_O2(i1:i2)) * 2.205;  % lbm/s

% dimensions of nozzle
D_center = 0.5; % inches
D_outer_2 = 0.30; % inches
L_2 = 5.375; % inches
A_port_2 = pi * D_center^2 / 4 + 8 * pi * D_outer_2^2 / 4; % cross_sectional area

% calculate r_dot and Go for fitting
delta_V_2 = mfuel_2 / rho; %change in Volume
delta_r_2 = sqrt(delta_V_2 / (pi * L_2));
r_dot_2 = delta_r_2 / t_2;   % in/sec

Go_2 = m_dot_O2_2 / A_port_2; % lbm/(s*in^2)

%% data point 3: 3rd Fire
load('bradyfire3.mat');
load('bradyfire2.mat');
i1 = start_index;
i2 = final_index;
t_3 = time(i2) - time(i1); % secs
mfuel_3 = mfuel / 10^3 * 2.205; % lbm;
m_dot_O2_3 = mean(m_dot_O2(i1:i2)) * 2.205;  % lbm/s

% dimensions of nozzle
D_center = 0.5; % inches
D_outer_3 = 0.30; % inches
D_outer_3_2 = 0.16;
L_3 = 5.375; % inches
A_port_3 = pi * D_center^2 / 4 + 7 * pi * D_outer_3^1 / 4 + 7 * pi * D_outer_3_2^2; % cross_sectional area

% calculate r_dot and Go for fitting
delta_V_3 = mfuel_3 / rho; %change in Volume
delta_r_3 = sqrt(delta_V_3 / (pi * L_3));
r_dot_3 = delta_r_3 / t_3;   % in/sec

Go_3 = m_dot_O2_3 / A_port_3; % lbm/(s*in^2)

%% Fit regression line
x = [Go_1 Go_2 Go_3 0.2385];
%Note: This is using data from "I LOVE DARREN" right now%
% x = [Go_1 Go_2 0.1844 0.2385];
% y = [r_dot_1 r_dot_2 0.1380 0.2076];
y = [r_dot_1 r_dot_2 r_dot_3 0.2076];

p = polyfit(log(x), log(y), 1);
n = p(1)
a = exp(p(2))

%% Predict mix ratio (using fire 2) 
L = 5.375;

%First Mix Ratio
p_1 = pi * D_center + 8 * pi * D_outer_1;
m_fuel_dot_1 = rho * L * p_1 * a * (m_dot_O2_1 / A_port_1) ^ n;
mixRatio_1 = m_dot_O2_1 / m_fuel_dot_1 * 10

%Second Mix ratio
p_2 = pi * D_center + 8 * pi * D_outer_2;
m_fuel_dot_2 = rho * L * p_2 * a * (m_dot_O2_2 / A_port_2) ^ n;
mixRatio_2 = m_dot_O2_2 / m_fuel_dot_2 * 10

%3rd Mix Ratio
p_3 = pi * D_center + 7 * pi * D_outer_3 + 7 * pi * D_outer_3_2;
m_fuel_dot_3 = rho * L * p_3 * a * (m_dot_O2_3 / A_port_3) ^ n;
mixRatio_3 = m_dot_O2_3 / m_fuel_dot_3 * 10

%%get surface areas
% SA_1 = p_1 * L
SA_2 = p_2 * L
SA_3 = p_3 * L
