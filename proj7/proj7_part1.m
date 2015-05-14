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
    % Combustor
    x1       = zeros(nsp,1);
    x1(iC2H4) = 1;  % Need to covert to molar ratio
    x1(iO2)  = phi(i);
    
    set(gas, 'T', Tref, 'P', P1, 'X', x1); 
    equilibrate(gas, 'HP');
    
    h1 = enthalpy_mass(gas);
    h_correct = h1 + (hf.CH2 - hf.C2H4) * massFraction(gas, iC2H4); %in a mass basis
    setState_HP(gas, [h_correct, P1]);
    equilibrate(gas, 'HP');
    
    %Combustor Outlet/ Nozzle Inlet
    To(i) = temperature(gas); 
    
    %% Frozen Flow
%     %Stagnation enthalpy is constant
%     h1 = enthalpy_mass(gas);
%     ho1 = h1; %assuming v1 about 0
%     
%     
    %Frozen Flow
    %Stagnation Enthalpy is constant
    % Stage One is before Nozzle
    % Stage 2 is at throat
    Ma = 1;
    P = P1;
    ho1 = h1;
    ho2 = 0;
    s1 = entropy_mass(gas); %%is this what we want?
    dT = 1;
    dP = 10;
    s2 = 0;
    
    
    while s2 < s1
        P = P - dP;
        ho2 = 0;
        T = To(i);
        while ho2 < ho1
            T = T - dT;
            if (T < 0)
                break;
            end
            set(gas, 'T', T, 'P', P);
            c = soundspeed(gas);
            h2 = enthalpy_mass(gas);
            V2 = Ma * c;
            ho2 = h2 + 0.5 * V2^2;
        end
        s2 = entropy_mass(gas)
    end
    
    
    T_t_frozen(i) = T;
    



    
    
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
figure;
plot(mixRatio, T_t_frozen);
xlabel('Mix Ratio');
ylabel('Nozzle Throat Temperature (K)');
title('Mix Ratio vs. Nozzle Throat Temperature');
set(gcf, 'color', 'white');

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

