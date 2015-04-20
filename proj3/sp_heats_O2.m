function [cp, cv, gamma, R] = sp_heats_O2(temp)

molarMass = 32; %g/mol
R = 8.314462 ./ molarMass * 10^3;

a = 25.48;
b = 1.520*10^-2;
c = -0.7155*10^-5;
d = 1.312*10^-9;

p = [d c b a];
cp = (polyval(p,temp)) * (1000 / molarMass);
cv = cp - R;
gamma = cp ./ cv;

end