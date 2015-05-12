clc
close all
clear

%Initialize gas object
gas = IdealGasMix('gri30.xml');

%Stagnation values
Po = pressure(gas);
To = temperature(gas);
a_o = soundspeed(gas); %speed of sound
rho_o = density(gas);

%Number of species in gas solution
nsp = nSpecies(gas);

%Find relevant species
iCH2 = speciesIndex(gas,'C2H4');
iO2  = speciesIndex(gas,'O2');
iCO2 = speciesIndex(gas,'CO2');
iH2O = speciesIndex(gas,'H2O');

%Declare initial mole fractions
x_r       = zeros(nsp,1);
x_r(iCH2) = 2;
x_r(iO2)  = 3;

x_p       = zeros(nsp,1);
x_p(iH2O) = 2;
x_p(iCO2) = 2;

%Initialize gas properties
set(gas,'T',To,'P',Po,'X',x_r);

%Equilibrate reaction
equilibrate(gas,'TP');

%Declare range of equivalence ratio
mixRatio = 1:0.1:10;

%Find heat of formation of CH2
h1 = enthalpy_mole(gas);

ao_over_ko = zeros(mixRatio,1);
for i = 1:length(mixRatio)
    R_hat = 0;
    Mo_hat = 0;
    ao_over_ko(i) = sqrt((R_hat * To) / (ko * Mo_hat));
end

%Calculate c*
c_star = (ao_over_ko) * (rho_o / rho_t) * (a_o / at);

