clc
close all

%molecular masses
MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2O = 18.016;

%enthalpy (J/mol) and entropy( J/(mol*K) ) of formation values
hf.H2O_vap = -241820;       sf.H2O_vap = 188.83;
hf.H2O_liq = -285830;       sf.H2O_liq = 69.92;
hf.O2 = 0;                  sf.O2 = 205.04;
hf.N2 = 0;                  sf.N2 = 191.61;
hf.H2 = 0;                  sf.H2 = 130.68;

%calculates the integrals
fun_O2_h = @(T)sp_heats(T,'O2');  %CHANGE sp_heats TO PER MOLS??? (I THINK IT IS PER KG RIGHT NOW)
fun_O2_s = @(T)sp_heats(T,'O2')./T;
fun_N2_h = @(T)sp_heats(T,'N2');
fun_N2_s = @(T)sp_heats(T,'N2')./T;
% fun_H2_h = @(T)sp_heats(T,'H2'); %DEFINE H2 IN SP_HEATS FUNCTION
% fun_H2_s = @(T)sp_heats(T,'H2')./T;%" " 
fun_H2O_vap_h = @(T)sp_heats(T,'H2O_vap'); %DEFINE IN SP_HEATS FUNCTION
fun_H2O_vap_s = @(T)sp_heats(T,'H2O_vap')./T;%" " 
% fun_H2O_liq_h = @(T)sp_heats(T,'H2O_liq'); %DEFINE IN SP_HEATS FUNCTION
% fun_H2O_liq_s = @(T)sp_heats(T,'H2O_liq')./T;%" " 

T_standard = 298;       P_standard = 1;
T = 300;                P = P_standard;
R = 8.314; %J/(mol*K) --->>????? is this the right "R"?????

% Calculate molar fractions for products [g/mol * g/mol] = [-]
lamda = 2;
phi = 1 ./ lamda;
N_prod.O2 = lamda * MM.O2;
N_prod.H2O = 2 * MM.H2O;
N_prod.N2 = lamda * 3.76 * MM.N2;
N_prod.sum = 2 * MM.H2O + lamda * 3.76 * MM.N2 + lamda * MM.O2;

y_prod.N2 = N_prod.N2 ./ N.sum;
y_prod.O2 = N_prod.O2 ./ N.sum;
y_prod.H2O = N_prodH2O ./ N.sum;


%calculate gibbs free energy of each species given balanced chemical
%reaction
g_prod.O2 = hf.O2 + integral(fun_O2_h, T_standard, T)...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T)) - R * log(P/P_standard));
g_prod.N2 =  hf.N2 + integral(fun_N2_h, T_standard, T)...
            - T * ((sf.N2 + integral(fun_N2_s, T_standard, T)) - R * log(P/P_standard));
% g_prod.H2 =  hf.H2 + integral(fun_H2_h, T_standard, T)...
%             - T * ((sf.H2 + integral(fun_H2_s, T_standard, T)) - R * log(P/P_standard))
g_prod.H2O_vap =  hf.H2O_vap + integral(fun_H2O_vap_h, T_standard, T)...
            - T * ((sf.H2O_vap + integral(fun_H2O_vap_s, T_standard, T)) - R * log(P/P_standard));

% Products:

g_prod.total = y_prod.N2 .* g_prod.N2 + y_prod.H2O .* g_prod.H2O + y_prod.O2 .* g_prod.O2;




%test
