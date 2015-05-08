function [h] = lucio_wgs_n(N_CO, N_H2O, N_CO2, N_H2, T)
% used to find enthalpys of all components in WGS reaction, 
% INPUT: mol fractions 
% OUTPUT: enthalpy for each element
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

N_prod.CO2 = N_CO2;
N_prod.H2 = N_H2;

m_prod.CO2 = N_prod.CO2 * MM.CO2;
m_prod.H2 = N_prod.H2 * MM.H2;


N_react.CO = N_CO;
N_react.H2O = N_H2O;

m_react.CO = N_react.CO * MM.CO;
m_react.H2O = N_react.H2O * MM.H2O;

% Enthalpy and Gibbs free energy of prod and react - all in J / kg
% react
H.H2O = (hf.H2O + integral(fun_H2O_h, T_standard, T)) * m_react.H2O / 1000;

H.CO = (hf.CO + integral(fun_CO_h, T_standard, T)) * m_react.CO / 1000;

% prod
H.H2 = (hf.H2 + integral(fun_H2_h, T_standard, T)) * m_prod.H2 / 1000;

H.CO2 = (hf.CO2 + integral(fun_CO2_h, T_standard, T)) * m_prod.CO2 / 1000;

h_CO = H.CO;
h_H2O = H.H2O;
h_CO2 = H.CO2;
h_H2 = H.H2;

h = h_CO + h_H2O + h_CO2 + h_H2;

end

