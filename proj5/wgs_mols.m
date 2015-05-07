function [N_CO, N_H2O, N_CO2, N_H2] = wgs_mols(T,P,Kp_wgs)
% used to find enthalpys of all components in WGS reaction,
% INPUT: T,P
% OUTPUT: mole fractions (N_i - N_n)
% CO + H2O -> CO2 + H2

% Universal gas constant R
R = 8.3144621;
T_standard = 298; %K
P_s = P; %Pa


syms N_CO N_H2O N_CO2 N_H2 N_total
assume(N_CO >= 0);
assume(N_H2O >= 0);
assume(N_CO2 >= 0);
assume(N_H2 >= 0);
assume(N_total >= 0);
eqns(1) = N_total == N_CO + N_H2O + N_CO2 + N_H2;
eqns(2) = Kp_wgs(T) == ((N_CO2 * N_H2) / (N_CO * N_H2O));
eqns(3) = 1 == N_CO + N_CO2;
eqns(4) = 10 == 2 * N_H2O + 2 * N_H2;
eqns(5) = 3 == N_CO + N_H2O + 2 * N_CO2;
S = solve(eqns, 'Real', true);
N_CO = double(S.N_CO);
N_H2O = double(S.N_H2O);
N_CO2 = double(S.N_CO2);
N_H2 = double(S.N_H2);




