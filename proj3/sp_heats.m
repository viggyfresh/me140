function [ cp, cv, gamma, R ] = sp_heats( temp , type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(type,'CO2')
    
    molarMass = 44.01; %grams/mole
    R = 8.314462 ./ molarMass * 10^3;
    a = 22.26;
    b = 5.981*10^-2;
    c = -3.501*10^-5;
    d = 7.469*10^-9;
    
elseif strcmp(type,'H2O')
    
    molarMass = 18.02;
    R = 8.314462 ./ molarMass * 10^3;
    a = 32.24;
    b = 0.1923*10^-2;
    c = 1.055*10^-5;
    d = -3.595*10^-9;
    
    
elseif strcmp(type,'N2')
    
    molarMass = 28.02;
    R = 8.314462 ./ molarMass * 10^3;
    a = 28.9;
    b = -0.1571*10^-2;
    c = 0.8081 * 10^-5;
    d = -2.873*10^-9;
    
elseif strcmp(type,'O2')
    
    molarMass = 32; %g/mol
    R = 8.314462 ./ molarMass * 10^3;
    a = 25.48;
    b = 1.520*10^-2;
    c = -0.7155*10^-5;
    d = 1.312*10^-9;
    
elseif strcmp(type,'air')
    
    R = 286.9;
    molarMass = 28.97;
    a = 28.11;
    b = 0.1967*10^-2;
    c = 0.4802*10^-5;
    d = -1.966*10^-9;

end

p = [d c b a];
cp = (polyval(p,temp)) * (1000 / molarMass);
cv = cp - R;
gamma = cp ./ cv;

end

