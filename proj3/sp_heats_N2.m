function [cp, cv, gamma, R] = sp_heats_N2(temp)

molarMass = 28.02;
R = 8.314462 ./ molarMass * 10^3;

a = 28.9;
b = -0.1571*10^-2;
c = 0.8081 * 10^-5;
d = -2.873*10^-9;

p = [d c b a];
cp = (polyval(p,temp)) * (1000 / molarMass);
cv = cp - R;
gamma = cp ./ cv;

end