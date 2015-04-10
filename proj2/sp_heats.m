function [cp,cv,gamma,R] = sp_heats(temp)

R=286.9;
molarMass=28.97;

a=28.11;
b=0.1967*10^-2;
c=0.4802*10^-5;
d=-1.966*10^-9;

p=[d c b a];
cp=(polyval(p,temp))*(1000/molarMass);
cv=cp-R;
gamma=cp./cv;

end

