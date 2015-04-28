function [ alpha ] = john(T, lambda)

P_standard = 101.325 .* 10^3;
P_sat = exp(-1.2914e8 ./ T^3 + 8.2048e5 ./ T^2 - 6522.8 ./ T + 25.5887);
N_a = (0.5 .* (lambda - 1) + (0.5 .* lambda .* 3.76));

alpha = (P_sat / P_standard * N_a) / (1 - P_sat / P_standard);

end

