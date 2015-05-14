clc
close all
clear

%%Part 1
%Initialize gas object
gas = IdealGasMix('gri30.xml');

%Stagnation values???
P1 = 6800000; %Pa
Tref = 298;


a_o = soundspeed(gas); % speed of sound
rho_o = density(gas);

%Number of species in gas solution
nsp = nSpecies(gas);

%Find relevant species
iC2H4 = speciesIndex(gas,'C2H4');
iO2  = speciesIndex(gas,'O2');
iCO2 = speciesIndex(gas,'CO2');
iH2O = speciesIndex(gas,'H2O');

% Declare range of equivalence ratios
mixRatio = 1:10;

% Molar masses (kg / mol)
MM.C2H4 = 28.052 / 1000;
MM.CH2 = MM.C2H4 / 2;
MM.H2O = 18.02 / 1000;
MM.CO2 = 44.01 / 1000;
MM.O2 = 32 / 1000;

% Heat of formations
hf.C2H4 = 52280 / MM.C2H4; % J/kg
hf.CO2 = -393520 / MM.CO2; % J/kg
hf.H2O =  -285830 / MM.H2O; % J/kg

% Find heat of formation of CH2
HHV_CH2 = 46.5 * 10^6; % J/kg
hf.CH2 = (2 * MM.CO2 * hf.CO2 + 2 * MM.H2O * hf.H2O + HHV_CH2 * 2 * MM.CH2) / (2 * MM.CH2)
phi = mixRatio * MM.C2H4 / MM.O2;

%%Combustor and Nozzle
for i=1:length(mixRatio)
    % Declare mole ratios
    x1       = zeros(nsp,1);
    x1(iC2H4) = 1;  % Need to covert to molar ratio
    x1(iO2)  = phi(i);
    
    set(gas, 'T', Tref, 'P', P1, 'X', x1); 
    equilibrate(gas, 'HP');
    
    h1 = enthalpy_mass(gas)
    h_correct = h1 + (hf.CH2 - hf.C2H4) * massFraction(gas, iC2H4) %in a mass basis
    setState_HP(gas, [h_correct, P1]);
    equilibrate(gas, 'HP');
    
    To(i) = temperature(gas)
    
    %% Frozen Flow
%     %Stagnation enthalpy is constant
%     h1 = enthalpy_mass(gas);
%     ho1 = h1; %assuming v1 about 0
%     
%     
%   
%   
%     
%     while h2 < h1
%         T =  T + dT;
%         set(gas, 'T', T, 'P', P1);
%         h2 = enthalpy_mass(gas);
%     end
%     
%     To_frozen(i) = T
    
    
    
%    % THIS IS WRONG!!!!!!!
%     s1 = entropy_mass(gas);   
%     s2 = 0;
%     T = Tref;
%     dT = 1;
%     
%     %Frozen Flow
%     %Stagnation Entropy is constant
%     while s2 < s1
%         T =  T + dT;
%         set(gas, 'T', T, 'P', Po);
%         s2 = entropy_mass(gas);
%     end
    
    
    %% Chemical Equilbrium
     
end



% 
figure;
plot(mixRatio, To);
xlabel('Mix Ratio');
ylabel('Nozzle Stagnation Temperature (K)');
title('Mix Ratio vs. Nozzle Stagnation Temperature');
set(gcf, 'color', 'white');
% 
% figure;
% plot(mixRatio, T_nozzle_throat_frozen);
% xlabel('Mix Ratio');
% ylabel('Nozzle Throat Temperature (K)');
% title('Mix Ratio vs. Nozzle Throat Temperature');
% set(gcf, 'color', 'white');
% 
% 
% %
% %
% % ao_over_ko = zeros(length(mixRatio),1);
% % for i = 1:length(mixRatio)
% %     R_hat = 0;
% %     Mo_hat = 0;
% %     ao_over_ko(i) = sqrt((R_hat * To) / (ko * Mo_hat));
% % end
% %
% % %Calculate c*
% % c_star = (ao_over_ko) * (rho_o / rho_t) * (a_o / at);

