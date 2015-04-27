clc;
clear all;
close all;

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
T_standard = 298;       P_standard = 101.3 * 10^3;
T_series = 25:5:1000;
T_series = T_series + 273;
P = P_standard;
R = 8.314;
alpha = 0;

% Calculate molar fractions for products [g/mol * g/mol] = [-]
lambda = 2;

for i=1:length(T_series)
    T = T_series(i);
    P_sat = exp(-1.2914e8 / T^3 + 8.2048e5 / T^2 - 6522.8 / T + 25.5887);
    
    % Saturation computations
    y_max = P_sat / P;
    N_a = (0.5 * (lambda - 1) + (0.5 * lambda * 3.76));
    y_test = 1 / (1 + N_a);
    
    if y_test > y_max
        beta = (y_max * N_a) / (1 - y_max);
        gamma = 1 - beta;
    else
        y_actual = y_test;
        beta = 1;
        gamma = 0;
    end
    
    m_prod.H2O_vap = beta * MM.H2O;
    m_prod.H2O_liq = gamma * MM.H2O;
    m_prod.H2O = MM.H2O;
    m_prod.O2 = 0.5 * (lambda - 1) * MM.O2;
    m_prod.N2 = 0.5 * lambda * 3.76 * MM.N2;
    m_prod.sum = m_prod.H2O + m_prod.O2 + m_prod.N2;
    
    mf_prod.H2O_vap = m_prod.H2O_vap ./ m_prod.sum;
    mf_prod.H2O_liq = m_prod.H2O_liq ./ m_prod.sum;
    mf_prod.H2O = m_prod.H2O ./ m_prod.sum;
    mf_prod.N2 = m_prod.N2 ./ m_prod.sum;
    mf_prod.O2 = m_prod.O2 ./ m_prod.sum;
    
    m_react.H2 = 1 * MM.H2;
    m_react.O2 = 0.5 * lambda * MM.O2;
    m_react.N2 = 0.5 * lambda * 3.76 * MM.N2;
    m_react.H2O_vap = alpha * MM.H2O;
    m_react.sum = m_react.H2 + m_react.O2 + m_react.N2 + m_react.H2O_vap;
    
    mf_react.H2 = m_react.H2 ./ m_react.sum;
    mf_react.O2 = m_react.O2 ./ m_react.sum;
    mf_react.N2 = m_react.N2 ./ m_react.sum;
    mf_react.H2O_vap = m_react.H2O_vap ./ m_react.sum;
    
    % Enthalpy and Gibbs free energy - all in J / kg
    H.O2 = hf.O2 + integral(fun_O2_h, T_standard, T);
    g.O2 = H.O2...
        - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    H.N2 = hf.N2 + integral(fun_N2_h, T_standard, T);
    g.N2 = H.N2...
        - T * ((sf.N2 + integral(fun_N2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    H.H2 = hf.H2 + integral(fun_H2_h, T_standard, T);
    g.H2 = H.H2...
        - T * ((sf.H2 + integral(fun_H2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    H.O2 = hf.O2 + integral(fun_O2_h, T_standard, T);
    g.O2 = H.O2...
        - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
        - R * log(P/P_standard));
    
    H.H2O_vap = hf.H2O_vap + integral(fun_H2O_vap_h, T_standard, T);
    g.H2O_vap = H.H2O_vap...
        - T * ((sf.H2O_vap + integral(fun_H2O_vap_s, T_standard, T))...
        - R * log(P/P_standard));
    
    H.H2O_liq = hf.H2O_liq + integral(fun_H2O_liq_h, T_standard, T);
    g.H2O_liq = H.H2O_liq...
        - T * ((sf.H2O_liq + integral(fun_H2O_liq_s, T_standard, T))...
        - R * log(P/P_standard));
    
    % Reactants & Products:
    g.react = mf_react.N2 .* g.N2 + mf_react.O2 .* g.O2...
        + mf_react.H2 .* g.H2 + mf_react.H2O_vap .* g.H2O_vap;
    
    % Lower heating value - all gas; higher heating value - all liquid
    g.LHV = mf_prod.N2 .* g.N2 + mf_prod.H2O .* g.H2O_vap...
        + mf_prod.O2 .* g.O2;
    g.HHV = mf_prod.N2 .* g.N2 + mf_prod.H2O .* g.H2O_liq...
        + mf_prod.O2 .* g.O2;
    g.actual = mf_prod.N2 .* g.N2 +  mf_prod.O2 .* g.O2...
        + mf_prod.H2O_vap .* g.H2O_vap + mf_prod.H2O_liq .* g.H2O_liq;
    
    % Calculate actual enthalpy change for denominator
    H_prod = mf_prod.N2 * H.N2 + mf_prod.O2 * H.O2...
        + mf_prod.H2O_vap * H.H2O_vap + mf_prod.H2O_liq * H.H2O_liq;
    H_react = mf_react.N2 * H.N2 + mf_react.O2 * H.O2...
        + mf_react.H2 * H.H2 + mf_react.H2O_vap * H.H2O_vap;
    H_actual = H_prod - H_react;
    
    
    % Calculate efficiencies
    efficiency.LHV(i) = -m_react.sum * (g.actual - g.react)...
        ./ (m_react.H2 * LHV);
    efficiency.HHV(i) =  -m_react.sum * (g.actual - g.react)...
        ./ (m_react.H2 * HHV);
    efficiency.actual(i) = (g.actual - g.react) ./ (m_react.H2 * H_actual);
end

% Something (might be) wrong with these values

% Compute Carnot efficiency for all temperatures
for i=1:length(T_series)
    T = T_series(i);
    n_carnot(i) = (1 - (T_standard / T)) * 100; %percent
end

% Plotting for Part 1
figure;
plot(T_series, efficiency.LHV * 100, T_series, efficiency.HHV * 100,...
    T_series, efficiency.actual * 100, T_series, n_carnot);
xlabel('Temperature (K)');
ylabel('Efficiency (%)');
legend('LHV', 'HHV', '\DeltaH', '\eta_c');
title('First Law Efficiency vs. Temperature');
set(gcf, 'color', 'w');
plotfixer;

%% Part 2

T_values = [80,220,650,800] + 273; % K
lambda_range = 1:1:10;
P_range = 1:1:40; % Pa

Rvar.O2 = R/MM.O2 * 10^3;
Rvar.H2 = R/MM.H2 * 10^3;
Rvar.N2 = R/MM.N2 * 10^3;
Rvar.H2O = R/MM.H2O * 10^3;

% change pressures, hold lamda
lambda = 2;
for i=1:length(T_values)
    T = T_values(i);
    for j = 1:length(P_range)
        % Enthalpy and Gibbs free energy - all in J / kg
        P = P_range(j) * 101.3 * 10^3;
        H.O2 = hf.O2 + integral(fun_O2_h, T_standard, T);
        g.O2 = H.O2...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
            - Rvar.O2 * log(P/P_standard));
        
        H.N2 = hf.N2 + integral(fun_N2_h, T_standard, T);
        g.N2 = H.N2...
            - T * ((sf.N2 + integral(fun_N2_s, T_standard, T))...
            - Rvar.N2 * log(P/P_standard));
        
        H.H2 = hf.H2 + integral(fun_H2_h, T_standard, T);
        g.H2 = H.H2...
            - T * ((sf.H2 + integral(fun_H2_s, T_standard, T))...
            - Rvar.H2 * log(P/P_standard));
        
        H.O2 = hf.O2 + integral(fun_O2_h, T_standard, T);
        g.O2 = H.O2...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
            - Rvar.O2 * log(P/P_standard));
        
        H.H2O_vap = hf.H2O_vap + integral(fun_H2O_vap_h, T_standard, T);
        g.H2O_vap = H.H2O_vap...
            - T * ((sf.H2O_vap + integral(fun_H2O_vap_s, T_standard, T))...
            - Rvar.H2O * log(P/P_standard));
        
        H.H2O_liq = hf.H2O_liq + integral(fun_H2O_liq_h, T_standard, T);
        g.H2O_liq = H.H2O_liq...
            - T * ((sf.H2O_liq + integral(fun_H2O_liq_s, T_standard, T))...
            - Rvar.H2O * log(P/P_standard));
        
        % Reactants & Products:
        g.react = mf_react.N2 .* g.N2 + mf_react.O2 .* g.O2...
            + mf_react.H2 .* g.H2 + mf_react.H2O_vap .* g.H2O_vap;
        
        % Lower heating value - all gas; higher heating value - all liquid
        g.LHV = mf_prod.N2 .* g.N2 + mf_prod.H2O .* g.H2O_vap...
            + mf_prod.O2 .* g.O2;
        g.actual = mf_prod.N2 .* g.N2 +  mf_prod.O2 .* g.O2...
            + mf_prod.H2O_vap .* g.H2O_vap + mf_prod.H2O_liq .* g.H2O_liq;
        
        % Calculate actual enthalpy change for denominator
        H_prod = mf_prod.N2 * H.N2 + mf_prod.O2 * H.O2...
            + mf_prod.H2O_vap * H.H2O_vap + mf_prod.H2O_liq * H.H2O_liq;
        H_react = mf_react.N2 * H.N2 + mf_react.O2 * H.O2...
            + mf_react.H2 * H.H2 + mf_react.H2O_vap * H.H2O_vap;
        H_actual = H_prod - H_react;
        
        
        % Calculate efficiencies
        eta_p(i, j) = -m_react.sum * (g.actual - g.react)...
            ./ (m_react.H2 * LHV);
    end
end

% Plotting for Part 2a
figure;
plot(P_range, eta_p(1, :) * 100, P_range, eta_p(2, :) * 100,...
     P_range, eta_p(3, :) * 100, P_range, eta_p(4, :) * 100);
xlabel('Pressure (Pa)');
ylabel('Efficiency (%)');
legend('80^{\circ}C', '220^{\circ}C', '650^{\circ}C', '800^{\circ}C');
title('First Law Efficiency vs. Pressure');
set(gcf, 'color', 'w');
plotfixer;

% Change lamda, keep pressure constant

%%

for i=1:length(T_values)
    T = T_values(i);
    for j = 1:length(lambda_range)
        % Enthalpy and Gibbs free energy - all in J / kg
        P = P_standard;
        lambda = lambda_range(j)
        m_prod.H2O_vap = beta * MM.H2O;
        m_prod.H2O_liq = gamma * MM.H2O;
        m_prod.H2O = MM.H2O;
        m_prod.O2 = 0.5 * (lambda - 1) * MM.O2;
        m_prod.N2 = 0.5 * lambda * 3.76 * MM.N2;
        m_prod.sum = m_prod.H2O + m_prod.O2 + m_prod.N2;
        
        mf_prod.H2O_vap = m_prod.H2O_vap ./ m_prod.sum;
        mf_prod.H2O_liq = m_prod.H2O_liq ./ m_prod.sum;
        mf_prod.H2O = m_prod.H2O ./ m_prod.sum;
        mf_prod.N2 = m_prod.N2 ./ m_prod.sum;
        mf_prod.O2 = m_prod.O2 ./ m_prod.sum;
        
        m_react.H2 = 1 * MM.H2;
        m_react.O2 = 0.5 * lambda * MM.O2;
        m_react.N2 = 0.5 * lambda * 3.76 * MM.N2;
        m_react.H2O_vap = alpha * MM.H2O;
        m_react.sum = m_react.H2 + m_react.O2 + m_react.N2 + m_react.H2O_vap;
        
        mf_react.H2 = m_react.H2 ./ m_react.sum;
        mf_react.O2 = m_react.O2 ./ m_react.sum;
        mf_react.N2 = m_react.N2 ./ m_react.sum;
        mf_react.H2O_vap = m_react.H2O_vap ./ m_react.sum;
        
        mf_prod
        
        H.O2 = hf.O2 + integral(fun_O2_h, T_standard, T);
        g.O2 = H.O2...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
            - Rvar.O2 * log(P/P_standard));
        
        H.N2 = hf.N2 + integral(fun_N2_h, T_standard, T);
        g.N2 = H.N2...
            - T * ((sf.N2 + integral(fun_N2_s, T_standard, T))...
            - Rvar.N2 * log(P/P_standard));
        
        H.H2 = hf.H2 + integral(fun_H2_h, T_standard, T);
        g.H2 = H.H2...
            - T * ((sf.H2 + integral(fun_H2_s, T_standard, T))...
            - Rvar.H2 * log(P/P_standard));
        
        H.O2 = hf.O2 + integral(fun_O2_h, T_standard, T);
        g.O2 = H.O2...
            - T * ((sf.O2 + integral(fun_O2_s, T_standard, T))...
            - Rvar.O2 * log(P/P_standard));
        
        H.H2O_vap = hf.H2O_vap + integral(fun_H2O_vap_h, T_standard, T);
        g.H2O_vap = H.H2O_vap...
            - T * ((sf.H2O_vap + integral(fun_H2O_vap_s, T_standard, T))...
            - Rvar.H2O * log(P/P_standard));
        
        H.H2O_liq = hf.H2O_liq + integral(fun_H2O_liq_h, T_standard, T);
        g.H2O_liq = H.H2O_liq...
            - T * ((sf.H2O_liq + integral(fun_H2O_liq_s, T_standard, T))...
            - Rvar.H2O * log(P/P_standard));
        
        % Reactants & Products:
        g.react = mf_react.N2 .* g.N2 + mf_react.O2 .* g.O2...
            + mf_react.H2 .* g.H2 + mf_react.H2O_vap .* g.H2O_vap;
        
        % Lower heating value - all gas; higher heating value - all liquid
        g.LHV = mf_prod.N2 .* g.N2 + mf_prod.H2O .* g.H2O_vap...
            + mf_prod.O2 .* g.O2;
        g.actual = mf_prod.N2 .* g.N2 +  mf_prod.O2 .* g.O2...
            + mf_prod.H2O_vap .* g.H2O_vap + mf_prod.H2O_liq .* g.H2O_liq;
        
        % Calculate actual enthalpy change for denominator
        H_prod = mf_prod.N2 * H.N2 + mf_prod.O2 * H.O2...
            + mf_prod.H2O_vap * H.H2O_vap + mf_prod.H2O_liq * H.H2O_liq;
        H_react = mf_react.N2 * H.N2 + mf_react.O2 * H.O2...
            + mf_react.H2 * H.H2 + mf_react.H2O_vap * H.H2O_vap;
        H_actual = H_prod - H_react;
        
        % Calculate efficiencies
        eta_lambda(i, j) = -m_react.sum * (g.actual - g.react)...
            ./ (m_react.H2 * LHV);
    end
end
figure
plot(lambda_range, eta_lambda(1, :) * 100, lambda_range, ...
    eta_lambda(2, :) * 100, lambda_range, eta_lambda(3, :)...
    * 100, lambda_range, eta_lambda(4, :) * 100)
