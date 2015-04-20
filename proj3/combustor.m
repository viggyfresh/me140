function [ T ] = combustor(MM, phi, To3)

To = 25 + 273.15; % Reference temperature
T_guess = To; % Initial guess

dT = 0.1; % Temperature step
lhv = 42800* MM.LHV;
target = phi*lhv + 17.85 * (deltaH_var_cp( To, To3, 'O2')*MM.O2/1000 + (79/21)...
    *deltaH_var_cp(To, To3, 'N2')*MM.N2/1000); %J/mol
right = 0;
coeff = (2 * 12.3 + 11.1) / 2;

while target > right
    T_guess = T_guess + dT;
    deltaH_CO2 =12.3 *phi * sp_heats_CO2(T_guess) * MM.CO2 / 1000* dT; 
    deltaH_H2O = 11.1*phi* sp_heats_H2O(T_guess)* MM.H2O / 1000* dT;
    deltaH_N2 = coeff * sp_heats_N2(T_guess)* MM.N2 / 1000 * dT;
    deltaH_O2 = 17.85 * (1-phi) * sp_heats_O2(T_guess)* MM.O2 / 1000 * dT;
    right = right + deltaH_CO2 + deltaH_H2O + deltaH_N2 + deltaH_O2;
    
    
end
T = T_guess;
end