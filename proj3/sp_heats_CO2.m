function [cp, cv, gamma, R] = sp_heats_CO2(temp)

molarMass = 44.01; %grams/mole
R = 8.314462 ./ molarMass * 10^3;

a = 22.26;
b = 5.981*10^-2;
c = -3.501*10^-5;
d = 7.469*10^-9;

p = [d c b a];
cp = (polyval(p,temp)) * (1000 / molarMass);
cv = cp - R;
gamma = cp ./ cv;

end