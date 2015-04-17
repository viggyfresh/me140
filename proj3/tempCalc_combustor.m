function [ T ] = tempCalc_combustor( hf )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

To = 25 + 273.15; % Reference temperature
T_guess = To; % Initial guess
dT = 0.001; % Temperature step
target=hf.JetA;
right = 12.3 * hf.CO2 + 11.1 * hf.H2O;

while abs(target-right)/target > 0.01
    T_guess = T_guess + dT;
    deltaH_CO2 = sp_heats_CO2(T_guess)/1000 * dT;
    deltaH_H2O = sp_heats_H2O(T_guess)/1000 * dT;
    deltaH_N2 = sp_heats_N2(T_guess)/1000 * dT;
    right = right + 12.3 * deltaH_CO2 + 11.1 * deltaH_H2O + 17.85 * (79/21) *...
        deltaH_N2;
    abs(target-right)/target
end
T = T_guess;
end

