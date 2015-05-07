function [N_CO, N_H2O, N_CO2, N_H2, total] = wgs_mols(P,Kp_wgs)
% used to find enthalpys of all components in WGS reaction,
% INPUT: T,P
% OUTPUT: mole fractions (N_i - N_n)
% CO + H2O -> CO2 + H2

% Universal gas constant R
R = 8.3144621;
T_s = 298; %K
P_s = P; %Pa


syms N_CO N_H2O N_CO2 N_H2 N_total
assume(N_CO >= 0);
assume(N_H2O >= 0);
assume(N_CO2 >= 0);
assume(N_H2 >= 0);
assume(N_total >= 0);
eqns(1) = N_total == N_CO + N_H2O + N_CO2 + N_H2;
eqns(2) = Kp_wgs == ((N_CO2 * N_H2) / (N_CO * N_H2O));
eqns(3) = 1 == N_CO + N_CO2;
eqns(4) = 10 == 2 * N_H2O + 2 * N_H2;
eqns(5) = 3 == N_CO + N_H2O + 2 * N_CO2;
S = solve(eqns, 'Real', true);
total = min(double(S.N_total));
N_CO = min(double(S.N_CO)) / total;
N_H2O = min(double(S.N_H2O)) / total;
N_CO2 = min(double(S.N_CO2)) / total;
N_H2 = min(double(S.N_H2)) / total;
