clc;
close all;
clear;

%% Part 1
% Initialize gas object
gas = IdealGasMix('me140_species.xml');

% Some values
P1 = 6800000; %Pa
Tref = 298;

% Declare range of mixture ratios
mixRatio = 1:0.25:10;

% Declare species indices
iC2H4 = speciesIndex(gas,'C2H4');
iO2  = speciesIndex(gas,'O2');
iCO2 = speciesIndex(gas,'CO2');
iH2O = speciesIndex(gas,'H2O');

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

%% Combustor and Nozzle
for i=1:length(mixRatio)
    [To(i), T_t_frozen(i)] = black_magic(gas, P1, phi(i), hf, 'frozen');
    [~, T_t_dissoc(i)] = black_magic(gas, P1, phi(i), hf, 'dissoc');
    X(:, i) = moleFractions(gas);
end

%% Plots
figure;
plot(mixRatio, To, mixRatio, T_t_frozen, mixRatio, T_t_dissoc);
xlabel('Mix Ratio');
ylabel('Temperature (K)');
title('Mix Ratio vs. Various Temperatures');
legend('T_0', 'T_t frozen', 'T_t dissociative');
set(gcf, 'color', 'white');
plotfixer;

figure;
plot(mixRatio, X(iC2H4,:),'g')
hold on;
plot(mixRatio, X(iO2,:),'b')
plot(mixRatio, X(iCO2,:),'r')
plot(mixRatio, X(iH2O,:),'c')
xlabel('Mix Ratio');
ylabel('Mole Fraction');
title('Mix Ratio vs. Mole Fractions');
legend('C_2H_4', 'O_2', 'CO_2', 'H_2O');
set(gcf, 'color', 'white');
plotfixer;

% ao_over_ko = zeros(length(mixRatio),1);
% for i = 1:length(mixRatio)
%     R_hat = 0;
%     Mo_hat = 0;
%     ao_over_ko(i) = sqrt((R_hat * To) / (ko * Mo_hat));
% end
%
% %Calculate c*
% c_star = (ao_over_ko) * (rho_o / rho_t) * (a_o / at);