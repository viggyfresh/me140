function [ T ] = tempCalc_combustor( hf )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

To = 25 + 273.15; % Reference temperature
T_guess = To; % Initial guess
dT = 0.0001; % Temperature step
target=hf.JetA;
right = 12.3 * hf.CO2 + 11.1 * hf.H2O; 

while abs((target-right)/target) > 0.1 %this tolerance value has to be high(weak), and dT low for loop to converge
    T_guess = T_guess + dT;
    target = target + 17.85 * sp_heats_air(T_guess)/1000 * (28.97/1000) * dT; %the (x/1000) converts it to a per molar basis
    deltaH_CO2 = 12.3 * sp_heats_CO2(T_guess)/1000 * (44.01/1000) * dT; 
    deltaH_H2O = 11.1 * sp_heats_H2O(T_guess)/1000 * (18.02/1000) * dT;
    deltaH_N2 = 17.85 * (79/21) * sp_heats_N2(T_guess)/1000 * (28.02/1000) * dT;
    right = right + 12.3 * deltaH_CO2 + 11.1 * deltaH_H2O + 17.85 * (79/21) *...
        deltaH_N2;
    abs((target-right)/target)
end
T = T_guess;
end

