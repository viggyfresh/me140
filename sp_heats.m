function [cp,cv,gamma,R] = sp_heats(temp)
%sp_heats.m takes in a range of temperatures and returns the thermodynamic
%properties cp, cv, and their ratio, gamma, in SI units
%   'temp': vector of temperatures. e.g. temp = 273:1800

%For air at STP
R=287;
molarMass=28.97;

%cp fit coefficients
a=28.11;
b=0.1967*10^-2;
c=0.4802*10^-5;
d=-1.966*10^-9;

%Use polyval to calculate cp. Divide by molar mass to obtain units of kJ/(kg*K)
p=[d c b a];
cp=(polyval(p,temp))*(1000/molarMass);
cv=cp-R;
gamma=cp./cv;

end

