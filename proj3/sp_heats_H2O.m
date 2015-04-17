function [cp, cv, gamma, R] = sp_heats_H2O(temp)

%This is for water vapor

molarMass = 18.02;
R = 8.314462 ./ molarMass * 10^3;

a = 32.24;
b = 0.1923*10^-2;
c = 1.055*10^-5;
d = -3.595*10^-9;

p = [d c b a];
cp = (polyval(p,temp)) * (1000 / molarMass);
cv = cp - R;
gamma = cp ./ cv;

end