function [ T ] = tempCalc_combustor( hf )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

To = 25 + 273.15; % Reference temperature
T_guess = To; % Initial guess
dT = 0.001; % Temperature step
target=hf.JetA;
right = 12.3 * hf.CO2 + 11.1 * hf.H2O; 

while abs((target-right)/target) > 0.1
    T_guess = T_guess + dT;
    % need to fix this line
    target = target + 17.85 * sp_heats_air(T_guess)/1000 * (28.97/1000) * dT; %the (x/1000) converts it to a per molar basis
    deltaH_CO2 = 12.3 * sp_heats_CO2(T_guess) * (44.01/1000) * dT; 
    deltaH_H2O = 11.1 * sp_heats_H2O(T_guess) * (18.02/1000) * dT;
    deltaH_N2 = 17.85 * (79/21) * sp_heats_N2(T_guess) * (28.02/1000) * dT;
    right = right + deltaH_CO2 + deltaH_H2O + deltaH_N2;
end
T = T_guess;
end

