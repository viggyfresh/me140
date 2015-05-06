function [deltaG_rxn,Kp] = lucio_reform(T)
% used to find deltaG_rxn of methane reformation process
% CH4 + H20 -> CO + 3*H2

% pressure of hydrogen fuel line
P_standard = 101325; %Pa
P = P_standard;

% Universal gas constant R
R = 8.3144621;

% Molecular masses - g/mol
MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2 = 2 * MM.H;
MM.H2O = 18.016;
MM.CH4 = MM.C + 4*MM.H;
MM.CO = MM.C + .5*MM.O2;
MM.air = 28.97;


% Enthalpy (J/kg) and entropy (J/(kg*K)) of formation values
% hf.H2O_vap = -241820 / MM.H2O * 1000;
% sf.H2O_vap = 188.83 / MM.H2O * 1000;
% hf.H2O_liq = -285830 / MM.H2O * 1000;
% sf.H2O_liq = 69.92 / MM.H2O * 1000;
hf.H2O = -285830 / MM.H2O * 1000; % based off liquid H2O
sf.H2O = 69.92 / MM.H2O * 1000; % based off liquid H2O
hf.O2 = 0;
sf.O2 = 205.04 / MM.O2 * 1000;
hf.N2 = 0;
sf.N2 = 191.61 / MM.N2 * 1000;
hf.H2 = 0;
sf.H2 = 130.68 / MM.H2 * 1000;
% added values 
hf.CH4 = -74873 / MM.CH4 * 1000; 
sf.CH4 = 186.251 / MM.CH4 * 1000;
hf.CO = -110527 / MM.CO * 1000;
sf.CO = 197.653 / MM.CO * 1000;

% Calculates the integrals
fun_H2_h = @(T)sp_heats(T,'H2');
fun_H2_s = @(T)sp_heats(T,'H2')./T;
fun_H2O_h = @(T)sp_heats(T,'H2O_vap');
fun_H2O_s = @(T)sp_heats(T,'H2O_vap')./T;
% added values
fun_CH4_h = @(T)sp_heats(T,'CH4');
fun_CH4_s = @(T)sp_heats(T,'CH4')./T;
fun_CO_h = @(T)sp_heats(T,'CO');
fun_CO_s = @(T)sp_heats(T,'CO')./T;

% Reference conditions
T_standard = 298;
P_standard = P;
P_sat = exp(-1.2914e8 / T^3 + 8.2048e5 / T^2 - 6522.8 / T + 25.5887);

% mols 
N_prod.H2 = 3;
N_prod.CO = 1;
N_prod.sum = N_prod.H2 + N_prod.CO;

y_prod.H2 = N_prod.H2 ./ N_prod.sum;
y_prod.CO = N_prod.CO ./ N_prod.sum;

m_prod.H2 = N_prod.H2 * MM.H2;
m_prod.CO = N_prod.CO * MM.CO;
m_prod.sum = m_prod.H2 + m_prod.CO;

mf_prod.H2 = m_prod.H2 ./ m_prod.sum;
mf_prod.CO = m_prod.CO ./ m_prod.sum;

N_react.H2O = 1;
N_react.CH4 = 1;
N_react.sum = N_react.H2O + N_react.CH4;

y_react.H2O = N_react.H2O ./ N_react.sum;
y_react.CH4 = N_react.CH4 ./ N_react.sum;

m_react.H2O = N_react.H2O * MM.H2O;
m_react.CH4 = N_react.CH4 * MM.CH4;
m_react.sum = m_react.H2O + m_react.CH4;

mf_react.H2O = m_react.H2O ./ m_react.sum;
mf_react.CH4 = m_react.CH4 ./ m_react.sum;

% Enthalpy and Gibbs free energy of prod and react - all in J / kg
H.H2O = hf.H2O + integral(fun_H2O_h, T_standard, T);
s_react.H2O = (sf.H2O + integral(fun_H2O_s, T_standard, T));
g_react.H2O = H.H2O - T * s_react.H2O;

H.CH4 = hf.CH4 + integral(fun_CH4_h, T_standard, T);
s_react.CH4 = (sf.CH4 + integral(fun_CH4_s, T_standard, T));
g_react.CH4 = H.CH4 - T * s_react.CH4;

H.H2 = hf.H2 + integral(fun_H2_h, T_standard, T);
s_prod.H2 = (sf.H2 + integral(fun_H2_s, T_standard, T));
g_prod.H2 = H.H2 - T * s_prod.H2;

H.CO = hf.CO + integral(fun_CO_h, T_standard, T);
s_prod.CO = (sf.CO + integral(fun_CO_s, T_standard, T));
g_prod.CO = H.CO - T * s_prod.CO;

% Reactants & Products:
g.react = (mf_react.CH4 .* g_react.CH4) + (mf_react.H2O .* g_react.H2O);

g.prod = (mf_prod.CO .* g_prod.CO) +  (mf_prod.H2 .* g_prod.H2);

% S.react = m_react.N2 .* s_react.N2 + m_react.O2 .* s_react.O2...
%     + m_react.H2 .* s_react.H2 + m_react.H2O_vap .* s_react.H2O_vap;
% S.react = S.react / 1000;
% S.prod = m_prod.N2 .* s_prod.N2 +  m_prod.O2 .* s_prod.O2...
%     + m_prod.H2O_vap .* s_prod.H2O_vap + m_prod.H2O_liq .* s_prod.H2O_liq;
% S.prod = S.prod / 1000;
% S_gen = S.prod - S.react;
% irrev = T * S_gen;

deltaG_rxn = (m_react.sum/1000) * (g.prod - g.react);

end

