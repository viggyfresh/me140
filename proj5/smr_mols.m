function [N_CO, N_H2O, N_CH4, N_H2] = smr_mols(T,P,Kp_smr)
% used to find enthalpys of all components in WGS reaction,
% INPUT: T,P
% OUTPUT: mole fractions (N_i - N_n) 
% CO + H2O -> CO2 + H2

% Universal gas constant R
R = 8.3144621;
T_standard = 298; %K
P_s = 101325; %Pa

% SOLVE system of equations
syms N_CO N_H2 N_CH4 N_H2O N_total
assume(N_CO >= 0);
assume(N_H2 >= 0);
assume(N_CH4 >= 0);
assume(N_H2O >= 0);
assume(N_total >= 0);
eqns(1) = N_total == N_CO + N_H2 + N_CH4 + N_H2O;
eqns(2) = Kp_smr == ((N_CO * N_H2^3) / (N_CH4 * N_H2O))...
    * ((P / P_s) / N_total)^2;
eqns(3) = 1 == N_CH4 + N_CO;
eqns(4) = 10 == 4 * N_CH4 + 2 * N_H2O + 2 * N_H2;
eqns(5) = 3 == N_H2O + N_CO;
S = solve(eqns, 'Real', true);
N_CO = double(S.N_CO);
N_H2 = double(S.N_H2);
N_CH4 = double(S.N_CH4);
N_H2O = double(S.N_H2O);





