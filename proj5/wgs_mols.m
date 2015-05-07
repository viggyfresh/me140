function [N_CO, N_H2O, N_CO2, N_H2] = wgs_mols(T,P)
% used to find enthalpys of all components in WGS reaction,
% INPUT: T,P
% OUTPUT: mole fractions (N_i - N_n) for all components
% CO + H2O -> CO2 + H2

% Universal gas constant R
R = 8.3144621;
T_standard = 298; %K
P_s = 101325; %Pa

% Molecular masses - g/mol
MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2 = 2 * MM.H;
MM.H2O = 18.016;
MM.CO2 = MM.C + MM.O2;
MM.CO = MM.C + .5*MM.O2;


syms N_CO positive N_H2 positive N_CH4 positive N_H2O positive N_total positive
eqns(1) = N_total == N_CO + N_H2 + N_CH4 + N_H2O;
eqns(2) = Kp_smr(j) == ((N_CO * N_H2^3) / (N_CH4 * N_H2O)) * ((P / P_s) / N_total)^2;
eqns(3) = 1 == N_CH4 + N_CO;
eqns(4) = 10 == 4 * N_CH4 + 2 * N_H2O + 2 * N_H2;
eqns(5) = 3 == N_H2O + N_CO;
S = solve(eqns, 'Real', true);
CO = double(S.N_CO);
H2 = double(S.N_H2);
CH4 = double(S.N_CH4);
H2O = double(S.N_H2O);
total = min(double(S.N_total));
B2.CO(i, j) = min(CO) / total;
B2.H2(i, j) = min(H2) / total;
B2.CH4(i, j) = min(CH4) / total;
B2.H2O(i, j) = min(H2O) / total;




