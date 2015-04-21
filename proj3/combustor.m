function [ T ] = combustor(MM, phi, To3)

To = 25 + 273.15; % Reference temperature
T_guess = To; % Initial guess

dT = 0.1; % Temperature step
lhv = 42800 * MM.JetA;
target = phi * lhv + 17.85 * (deltaH_var_cp(To, To3, 'O2', 1, phi, MM) * MM.O2 / 1000 ...
         + (79/21) * deltaH_var_cp(To, To3, 'N2', 1, phi, MM) * MM.N2 / 1000); %J/mol
right = 0;
coeff = (2 * 12.3 + 11.1) / 2;

while target > right
    T_guess = T_guess + dT;
    deltaH_CO2 = 12.3 * phi * sp_heats(T_guess, 'CO2') * MM.CO2 / 1000 * dT; 
    deltaH_H2O = 11.1 * phi * sp_heats(T_guess, 'H2O') * MM.H2O / 1000 * dT;
    deltaH_N2 = coeff * (79 / 21) * sp_heats(T_guess, 'N2') * MM.N2 / 1000 * dT;
    deltaH_O2 = coeff * (1 - phi) * sp_heats(T_guess, 'O2') * MM.O2 / 1000 * dT;
    right = right + deltaH_CO2 + deltaH_H2O + deltaH_N2 + deltaH_O2;
end
T = T_guess;
end