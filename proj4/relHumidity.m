function [ alpha, RH ] = relHumidity( T, lambda )

P_standard = 101.3 .* 10^3;
P_sat = exp(-1.2914e8 ./ T^3 + 8.2048e5 ./ T^2 - 6522.8 ./ T + 25.5887);

N_a = (0.5 .* (lambda - 1) + (0.5 .* lambda .* 3.76));
y_max = P_sat / P_standard;
y_test = 0;
alpha = 0;
dalpha = .1;
N_tot_H2O = 0;

while y_test < y_max
    N_tot_H2O = 1+alpha; % 1+ alpha = beta + gamma
    y_test = N_tot_H2O ./ (N_tot_H2O + N_a);
    alpha = alpha + dalpha
end


RH = ((1+alpha)./(N_tot_H2O + N_a)).* P_standard ./ P_sat;
%RH = P_v / P_sat_35C;
end

