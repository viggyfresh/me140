function [h] = lucio_smr_n(N_CH4, N_H2O, N_CO, N_H2, T)
% used to find enthalpies of reformation (SMR) reaction
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
hf.H2O = -241820 / MM.H2O * 1000;
sf.H2O = 188.83 / MM.H2O * 1000;
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

% mols 
N_prod.H2 = N_H2;
N_prod.CO = N_CO;

m_prod.H2 = N_prod.H2 * MM.H2;
m_prod.CO = N_prod.CO * MM.CO;

N_react.H2O = N_H2O;
N_react.CH4 = N_CH4;

m_react.H2O = N_react.H2O * MM.H2O;
m_react.CH4 = N_react.CH4 * MM.CH4;

% Enthalpy and Gibbs free energy of prod and react - all in J / kg
H.H2O = (hf.H2O + integral(fun_H2O_h, T_standard, T)) * m_react.H2O / 1000;

H.CH4 = (hf.CH4 + integral(fun_CH4_h, T_standard, T)) * m_react.CH4 / 1000;

H.H2 = (hf.H2 + integral(fun_H2_h, T_standard, T)) * m_prod.H2 / 1000;

H.CO = (hf.CO + integral(fun_CO_h, T_standard, T)) * m_prod.CO / 1000;

h_CH4 = H.CH4;
h_H2O = H.H2O;
h_CO = H.CO;
h_H2 = H.H2;

h = h_CH4 + h_H2O + h_CO + h_H2;
end

