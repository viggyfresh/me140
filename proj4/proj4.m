clc
close all

%molecular masses
MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2 = 2*MM.H;

MM.H2O = 18.016;

%enthalpy (J/mol) and entropy( J/(mol*K) ) of formation values
hf.H2O_vap = -241820;       sf.H2O_vap = 188.83;
hf.H2O_liq = -285830;       sf.H2O_liq = 69.92;
hf.O2 = 0;                  sf.O2 = 205.04;
hf.N2 = 0;                  sf.N2 = 191.61;
hf.H2 = 0;                  sf.H2 = 130.68;

%calculates the integrals
fun_O2_h = @(T)sp_heats(T,'O2');
fun_O2_s = @(T)sp_heats(T,'O2')./T;
fun_N2_h = @(T)sp_heats(T,'N2');
fun_N2_s = @(T)sp_heats(T,'N2')./T;
fun_H2_h = @(T)sp_heats(T,'H2');
fun_H2_s = @(T)sp_heats(T,'H2')./T;
fun_H2O_vap_h = @(T)sp_heats(T,'H2O_vap');
fun_H2O_vap_s = @(T)sp_heats(T,'H2O_vap')./T;
fun_H2O_liq_h = @(T)sp_heats(T,'H2O_liq');
fun_H2O_liq_s = @(T)sp_heats(T,'H2O_liq')./T;

T_standard = 298;       P_standard = 1;
T = 300;                P = P_standard;
R = 8.314; %J/(mol*K) --->>????? is this the right "R"?????

% Calculate molar fractions for products [g/mol * g/mol] = [-]
lambda = 2;

% Lower heating value - all gas
N_prod.H2O_vap = 1;
N_prod.H2O_liq = 0;
N_prod.O2 = 0.5 * (lambda - 1);
N_prod.N2 = 0.5 * lambda * 3.76;
N_prod.sum = N_prod.H2O_vap + N_prod.H2O_liq + N_prod.O2 + N_prod.N2;

y_prod.H2O_vap = N_prod.H2O_vap ./ N_prod.sum;
y_prod.H2O_liq = N_prod.H2O_liq ./ N_prod.sum;
y_prod.N2 = N_prod.N2 ./ N_prod.sum;
y_prod.O2 = N_prod.O2 ./ N_prod.sum;

N_react.H2 = 1;
N_react.O2 = 0.5 * lambda;
N_react.N2 = 0.5 * lambda * 3.76;
N_react.H2O_vap = 0;
N_react.sum = N_react.H2 + N_react.O2 + N_react.N2 + N_react.H2O_vap;

y_react.H2 = N_react.H2 ./ N_react.sum;
y_react.O2 = N_react.O2 ./ N_react.sum;
y_react.N2 = N_react.N2 ./ N_react.sum;
y_react.H2O_vap = N_react.H2O_vap ./ N_react.sum;

%calculate gibbs free energy of each species given balanced chemical
%reaction
% Gibbs free energy of reactants - all in J / kg
g_react.O2 = hf.O2 + integral(fun_O2_h, T_standard, T)...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T)) - R * log(P/P_standard));
        
g_react.N2 =  hf.N2 + integral(fun_N2_h, T_standard, T)...
            - T * ((sf.N2 + integral(fun_N2_s, T_standard, T)) - R * log(P/P_standard));
        
g_react.H2 =  hf.H2 + integral(fun_H2_h, T_standard, T)...
            - T * ((sf.H2 + integral(fun_H2_s, T_standard, T)) - R * log(P/P_standard));
        
% Gibbs free energy of products 
g_prod.O2 = hf.O2 + integral(fun_O2_h, T_standard, T)...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T)) - R * log(P/P_standard));
g_prod.N2 =  hf.N2 + integral(fun_N2_h, T_standard, T)...
            - T * ((sf.N2 + integral(fun_N2_s, T_standard, T)) - R * log(P/P_standard));
% g_prod.H2 =  hf.H2 + integral(fun_H2_h, T_standard, T)...
%             - T * ((sf.H2 + integral(fun_H2_s, T_standard, T)) - R * log(P/P_standard))
g_prod.H2O_vap =  hf.H2O_vap + integral(fun_H2O_vap_h, T_standard, T)...
            - T * ((sf.H2O_vap + integral(fun_H2O_vap_s, T_standard, T)) - R * log(P/P_standard));

% Reactants & Products:
g_react.total = y_react.N2 .* g_react.N2 + y_react.O2 .* g_react.O2 + y_react.H2 .* g_react.H2

g_prod.total = y_prod.N2 .* g_prod.N2 + y_prod.H2O .* g_prod.H2O_vap + y_prod.O2 .* g_prod.O2



%test
