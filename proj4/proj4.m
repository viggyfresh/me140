clc
clear all
close all

% Molecular masses - grams
MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2 = 2*MM.H;
MM.H2O = 18.016;
MM.CO2 = MM.C + MM.O2;
MM.air = 28.97; 

% Enthalpy (J/kg) and entropy (J/(kg*K)) of formation values
hf.H2O_vap = -241820 / MM.H2O * 1000;
sf.H2O_vap = 188.83 / MM.H2O * 1000;
hf.H2O_liq = -285830 / MM.H2O * 1000;
sf.H2O_liq = 69.92 / MM.H2O * 1000;
hf.O2 = 0;
sf.O2 = 205.04 / MM.O2 * 1000;
hf.N2 = 0;
sf.N2 = 191.61 / MM.N2 * 1000;
hf.H2 = 0;
sf.H2 = 130.68 / MM.H2 * 1000;

% Fuel heating values for H2 (J/kg)
LHV = 120.21 * 10^6;
HHV = 142.18 * 10^6;

% Calculates the integrals
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

% Reference conditions
T_standard = 298;       P_standard = 1;
T_series = 25:5:1000;
T_series = T_series + 273;
P = P_standard;
R = 8.314;

% Calculate molar fractions for products [g/mol * g/mol] = [-]
lambda = 2;

for i=1:length(T_series)
    T = T_series(i);
    
    % Lower heating value - all gas; higher heating value - all liquid
    m_prod.H2O_vap = 1 * MM.H2O;
    m_prod.H2O_liq = 0 * MM.H2O;
    m_prod.O2 = 0.5 * (lambda - 1) * MM.O2;
    m_prod.N2 = 0.5 * lambda * 3.76 * MM.N2;
    m_prod.sum = m_prod.H2O_vap + m_prod.H2O_liq + m_prod.O2 + m_prod.N2;
    
    mf_prod.H2O_vap = m_prod.H2O_vap ./ m_prod.sum;
    mf_prod.H2O_liq = m_prod.H2O_liq ./ m_prod.sum;
    mf_prod.N2 = m_prod.N2 ./ m_prod.sum;
    mf_prod.O2 = m_prod.O2 ./ m_prod.sum;
    
    m_react.H2 = 1 * MM.H2;
    m_react.O2 = 0.5 * lambda * MM.O2;
    m_react.N2 = 0.5 * lambda * 3.76 * MM.N2;
    m_react.H2O_vap = 0;
    m_react.sum = m_react.H2 + m_react.O2 + m_react.N2 + m_react.H2O_vap;
    
    mf_react.H2 = m_react.H2 ./ m_react.sum;
    mf_react.O2 = m_react.O2 ./ m_react.sum;
    mf_react.N2 = m_react.N2 ./ m_react.sum;
    mf_react.H2O_vap = m_react.H2O_vap ./ m_react.sum;
    
    % Gibbs free energy - all in J / kg
    g.O2 = hf.O2 + integral(fun_O2_h, T_standard, T)...
        - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    g.N2 =  hf.N2 + integral(fun_N2_h, T_standard, T)...
        - T * ((sf.N2 + integral(fun_N2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    g.H2 =  hf.H2 + integral(fun_H2_h, T_standard, T)...
        - T * ((sf.H2 + integral(fun_H2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    g.O2 = hf.O2 + integral(fun_O2_h, T_standard, T)...
        - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    g.N2 =  hf.N2 + integral(fun_N2_h, T_standard, T)...
        - T * ((sf.N2 + integral(fun_N2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    g.H2O_vap =  hf.H2O_vap + integral(fun_H2O_vap_h, T_standard, T)...
        - T * ((sf.H2O_vap + integral(fun_H2O_vap_s, T_standard, T))...
        - R * log(P/P_standard));
    
    g.H2O_liq =  hf.H2O_liq + integral(fun_H2O_liq_h, T_standard, T)...
        - T * ((sf.H2O_liq + integral(fun_H2O_liq_s, T_standard, T))...
        - R * log(P/P_standard));
    
    % Reactants & Products:
    g.total = mf_react.N2 .* g.N2 + mf_react.O2 .* g.O2...
        + mf_react.H2 .* g.H2;
    g.LHV = mf_prod.N2 .* g.N2 + mf_prod.H2O_vap .* g.H2O_vap...
        + mf_prod.O2 .* g.O2;
    % we use mf_prod.H2O_vap again because coefficient is the same
    g.HHV = mf_prod.N2 .* g.N2 + mf_prod.H2O_vap .* g.H2O_liq...
        + mf_prod.O2 .* g.O2;
    
    %calculate efficiencies
    efficiency.LHV(i) = -m_react.sum * (g.LHV - g.total)...
        ./ (m_react.H2 * LHV);
    efficiency.HHV(i) =  -m_react.sum * (g.HHV - g.total)...
        ./ (m_react.H2 * HHV);
end

% Something (might be) wrong with these values
figure;
plot(T_series, efficiency.LHV * 100, T_series, efficiency.HHV * 100);
xlabel('Temperature (K)');
ylabel('Efficiency (%)');
legend('LHV', 'HHV');
title('First Law Efficiency vs. Temperature');
set(gcf, 'color', 'w');
plotfixer;

%% first law efficiency (actual values)

%% carnot efficiency 

for i=1:length(T_series)
    T = T_series(i);
    n_carnot(i) = (1 - (T_standard / T))*100; %percent
end
hold on
plot(T_series,n_carnot,'g');


%% question 2 

T_values = [80,220,650,80]; % celsius
lambda_range = 1:.1:10;
P_range = 1:.1:40; %Pa

% change pressures, hold lamda
lamda = 2;


