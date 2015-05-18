clc;
clear all;
load('bradycheated.mat');
% i1 = 38;
% i2 = 6384;
i1 = 2403;
i2 = 5046;

Po = mean(chamP(i1:i2)) * 1000 + 101325;
D = 0.605; %inches
D = D / 39.370; %meters
At = pi * D^2 / 4;
t = time(i2) - time(i1) %secs
m_dot_fuel = mfuel / 10^3 / t; %kg
m_dot_O2 = mean(m_dot_O2(i1:i2)); 
mdot = m_dot_fuel + m_dot_O2;

cstar = Po / (mdot / At)
mixRatio = m_dot_O2 / m_dot_fuel