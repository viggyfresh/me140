function [deltaG_rxn] = lucio_wgs(T)
% used to find deltaG_rxn of water-gas-steam reaction
% CO + H2O -> CO2 + H2

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
MM.CO2 = MM.C + MM.O2;
MM.CO = MM.C + .5*MM.O2;
MM.air = 28.97;

Rvar.O2 = R / MM.O2 * 10^3;
Rvar.H2 = R / MM.H2 * 10^3;
Rvar.H2O = R / MM.H2O * 10^3;
Rvar.CO2 = R / MM.CO2 * 10^3;
Rvar.CO = R / MM.CO * 10^3;

% Enthalpy (J/kg) and entropy (J/(kg*K)) of formation values
% hf.H2O_vap = -241820 / MM.H2O * 1000;
% sf.H2O_vap = 188.83 / MM.H2O * 1000;
% hf.H2O_liq = -285830 / MM.H2O * 1000;
% sf.H2O_liq = 69.92 / MM.H2O * 1000;
hf.H2O = -241820 / MM.H2O * 1000;
sf.H2O = 188.83 / MM.H2O * 1000;
hf.H2 = 0;
sf.H2 = 130.68 / MM.H2 * 1000;
hf.CO = -110527 / MM.CO * 1000;
sf.CO = 197.653 / MM.CO * 1000;
hf.CO2 = -393522 / MM.CO2 * 1000;
sf.CO2 = 213.795 / MM.CO2 * 1000;


% Calculates the integrals
fun_H2_h = @(T)sp_heats(T,'H2');
fun_H2_s = @(T)sp_heats(T,'H2')./T;
fun_H2O_h = @(T)sp_heats(T,'H2O_vap'); % assume all vapor
fun_H2O_s = @(T)sp_heats(T,'H2O_vap')./T; % assue all vapor
fun_CO_h = @(T)sp_heats(T,'CO');
fun_CO_s = @(T)sp_heats(T,'CO')./T;
fun_CO2_h = @(T)sp_heats(T,'CO2');
fun_CO2_s = @(T)sp_heats(T,'CO2')./T;

% Reference conditions
T_standard = 298;
P_standard = P;
P_sat = exp(-1.2914e8 / T^3 + 8.2048e5 / T^2 - 6522.8 / T + 25.5887);

N_prod.CO2 = 1;
N_prod.H2 = 1;
N_prod.sum = N_prod.CO2 + N_prod.H2;

m_prod.CO2 = N_prod.CO2 * MM.CO2;
m_prod.H2 = N_prod.H2 * MM.H2;
m_prod.sum = m_prod.CO2 + m_prod.H2;

mf_prod.H2 = m_prod.H2 ./ m_prod.sum;
mf_prod.CO2 = m_prod.CO2 ./ m_prod.sum;

N_react.CO = 1;
N_react.H2O = 1;

m_react.CO = N_react.CO * MM.CO;
m_react.H2O = N_react.H2O * MM.H2O;
m_react.sum = m_react.CO + m_react.H2O;

mf_react.CO = m_react.CO ./ m_react.sum;
mf_react.H2O = m_react.H2O ./ m_react.sum;
mf_prod.H2 = m_prod.H2 ./ m_prod.sum;
mf_prod.CO2 = m_prod.CO2 ./ m_prod.sum;

% Enthalpy and Gibbs free energy of prod and react - all in J / kg
% react
H.H2O = hf.H2O + integral(fun_H2O_h, T_standard, T);
s_react.H2O = (sf.H2O + integral(fun_H2O_s, T_standard, T));
g_react.H2O = H.H2O - T * s_react.H2O;

H.CO = hf.CO + integral(fun_CO_h, T_standard, T);
s_react.CO = (sf.CO + integral(fun_CO_s, T_standard, T));
g_react.CO = H.CO - T * s_react.CO;

% prod
H.H2 = hf.H2 + integral(fun_H2_h, T_standard, T);
s_prod.H2 = (sf.H2 + integral(fun_H2_s, T_standard, T));
g_prod.H2 = H.H2 - T * s_prod.H2;

H.CO2 = hf.CO2 + integral(fun_CO2_h, T_standard, T);
s_prod.CO2 = (sf.CO2 + integral(fun_CO2_s, T_standard, T));
g_prod.CO2 = H.CO2 - T * s_prod.CO2;

% Reactants & Products:
g.react = (mf_react.H2O .* g_react.H2O) + (mf_react.CO .* g_react.CO);
g.prod = (mf_prod.H2 .* g_prod.H2) +  (mf_prod.CO2 .* g_prod.CO2);

deltaG_rxn = (m_react.sum/1000) * (g.prod - g.react);

end

