clc;
clear all;
load('bradyfire1.mat');
% i1 = 38;
% i2 = 6384;
i1 = start_index;
i2 = final_index;

plot(time, chamP, 'r', time, 10*thrust, 'b:', time, resP, 'k--', time, P_choked_min, 'g-')
legend('Chamber [kPa gage]', '10*Thrust [N]', 'Reservoir [kPa gage]', 'Unchoked orifice pressure')

Po = mean(chamP(i1:i2)) * 1000 + 101325
%D = 0.605; %inches, original
D = 23/64; %first fire
D = D / 39.370; %meters
At = pi * D^2 / 4;
t = time(i2) - time(i1) %secs
m_dot_fuel = mfuel / 10^3 / t; %kg
m_dot_O2 = mean(m_dot_O2(i1:i2)); 
mdot = m_dot_fuel + m_dot_O2;

cstar = Po / (mdot / At)
mixRatio = m_dot_O2 / m_dot_fuel