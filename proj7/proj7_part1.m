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

% %Initialize gas properties
% set(gas, 'T', To, 'P', Po, 'X', x_r);
%
% %Equilibrate reaction
% equilibrate(gas,'TP');

% Declare range of equivalence ratios
mixRatio = 1:0.1:10;

% Molar masses (kg / mol)
MM.C2H4 = 28.06 / 1000;
MM.CH2 = MM.C2H4 / 2;
MM.H2O = 18.02 / 1000;
MM.CO2 = 44.01 / 1000;

% Heat of formations
hf.C2H4 = 52280 / MM.C2H4; % J/kg
hf.CO2 = -393520 / MM.CO2; % J/kg
hf.H2O =  -285830 / MM.H2O; % J/kg

% Find heat of formation of CH2
HHV_CH2 = 46.5 * 10^6; % J/kg
hf.CH2 = (2 * MM.CO2 * hf.CO2 + 2 * MM.H2O * hf.H2O + HHV_CH2 * 2 * MM.CH2) / (2 * MM.CH2)
% set(gas, 'T', To, 'P', Po, 'X', x_r);
% equilibrate(gas, 'TP');

for i=1:length(mixRatio)
    
    h1 = enthalpy_mass(gas);
    h_correct = h1 + (hf.CH2 - hf.C2H4) * massFraction(gas, iCH2);
    setState_HP(gas, [h_correct, Po]);
    equilibrate(gas, 'HP');
    
    
    
    
end




% h_r = enthalpy_mole(gas);
%
% h_p = 0;
% T = To;
% dT = 0.1;
%
% while h_p < h_r
%     T = T + dT;
%     set(gas, 'T', T, 'P', P, 'X', x_p);
%     h_p = enthalpy_mole(gas);
% end
%
%
% ao_over_ko = zeros(length(mixRatio),1);
% for i = 1:length(mixRatio)
%     R_hat = 0;
%     Mo_hat = 0;
%     ao_over_ko(i) = sqrt((R_hat * To) / (ko * Mo_hat));
% end
%
% %Calculate c*
% c_star = (ao_over_ko) * (rho_o / rho_t) * (a_o / at);

