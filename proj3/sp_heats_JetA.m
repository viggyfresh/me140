function [cp, cv, gamma, R] = sp_heats_JetA(T, phi, MM)

[cp_CO2, ~, ~, R_CO2] = sp_heats(T, 'CO2'); % J / kg * K
[cp_H2O, ~, ~, R_H2O] = sp_heats(T, 'H2O'); % J / kg * K
[cp_N2, ~, ~, R_N2] = sp_heats(T, 'N2'); % J / kg * K
[cp_O2, ~, ~, R_O2] = sp_heats(T, 'O2'); % J / kg * K

coeff = (2 * 12.3 + 11.1) / 2;

%Mass quantities
M_total = phi .* (12.3 * MM.CO2 + 11.1 * MM.H2O) + ...
    coeff * ((79 / 21) * MM.N2 + (1 - phi) .* MM.O2); %g
m_CO2 = phi .* 12.3 * MM.CO2 ./ M_total;
m_H2O = phi .* 11.1 * MM.H2O ./ M_total;
m_N2 = coeff * (79 / 21) * MM.N2 ./ M_total;
m_O2 = coeff .* (1 - phi) .* MM.O2 ./ M_total;

cp = m_CO2 .* cp_CO2 + m_H2O .* cp_H2O + m_N2 .* cp_N2 + m_O2 .* cp_O2; % J / kg * K

%Lumped R
R = m_CO2 .* R_CO2 + m_H2O .* R_H2O + m_N2 .* R_N2 + m_O2 .* R_O2;

% Cv and gamma
cv = cp - R; % J / kg * K
gamma = cp ./ cv;


end

