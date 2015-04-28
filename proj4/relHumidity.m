function [ alpha, RH ] = relHumidity( T, lambda )

P_standard = 101.325 .* 10^3;
P_sat = exp(-1.2914e8 ./ T^3 + 8.2048e5 ./ T^2 - 6522.8 ./ T + 25.5887);

N_a = (0.5 .* (lambda - 1) + (0.5 .* lambda .* 3.76));
y_max = P_sat / P_standard;
alpha = 0;
dalpha = .01;
N_tot_H2O = 1 + alpha; %(fully saturated)
y_test = N_tot_H2O / (N_tot_H2O + N_a);

if y_test < y_max
    y_test = 0;
    while y_test < y_max
        N_tot_H2O = 1 + alpha; % 1 + alpha = beta + gamma
        y_test = N_tot_H2O ./ (N_tot_H2O + N_a);
        alpha = alpha + dalpha;
    end
    
else
    alpha = 0;
end

RH = ((alpha)./(alpha + N_a)).* P_standard ./ P_sat;

end

