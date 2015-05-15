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
mixRatio = 1:10;

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
tic
for i=1:length(mixRatio)
    i
    [To(i), T_t_frozen(i), c_star_frozen(i), T_e_frozen(i), V_e_frozen(i)] = black_magic(gas, P1, phi(i), hf, 'frozen');
    [~, T_t_dissoc(i), c_star_dissoc(i), T_e_dissoc(i), V_e_dissoc(i)] = black_magic(gas, P1, phi(i), hf, 'dissoc');
    X(:, i) = moleFractions(gas);
end
toc

%% Plots

% Plot of mole ratios
% figure;
% plot(mixRatio, X(iC2H4,:),'g')
% hold on;
% plot(mixRatio, X(iO2,:),'b')
% plot(mixRatio, X(iCO2,:),'r')
% plot(mixRatio, X(iH2O,:),'c')
% xlabel('Mixture Ratio');
% ylabel('Mole Fraction');
% title('Mixture Ratio vs. Mole Fractions');
% legend('C_2H_4', 'O_2', 'CO_2', 'H_2O');
% set(gcf, 'color', 'white');
% plotfixer;

% Plot of Throat temperature and stag temperature
figure;
plot(mixRatio, To, 'r', mixRatio, T_t_frozen, '--b', mixRatio, T_t_dissoc, 'b', mixRatio, T_e_frozen, '--g', mixRatio, T_e_dissoc, 'g');
xlabel('Mixture Ratio');
ylabel('Temperature (K)');
title('Mixture Ratio vs. Various Temperatures');
%legend('T_0', 'T_t frozen', 'T_t dissociative', 'T_e frozen', 'T_e dissociative');
set(gcf, 'color', 'white');
plotfixer;

% Plot of c star
hold on;
plot(mixRatio, c_star_frozen, '--k', mixRatio, c_star_dissoc, 'k');
xlabel('Mixture Ratio');
ylabel('c^* (m/s)');
title('c^* vs. Mixture Ratio');
% legend('Frozen', 'Dissociative');
ylim([0 6000])
set(gcf, 'color', 'white');

%Plot of Velocity
hold on;
plot(mixRatio, V_e_frozen, '--m', mixRatio, V_e_dissoc, 'm');

legend('T_0', 'T_t frozen', 'T_t dissociative', 'T_e frozen', 'T_e dissociative', ...
    'C^* frozen', 'c^* dissociative', 'V_e frozen', 'V_e dissociative');

plotfixer;

% ao_over_ko = zeros(length(mixRatio),1);
%     R_hat = 0;
%     Mo_hat = 0;
%     ao_over_ko(i) = sqrt((R_hat * To) / (ko * Mo_hat));
% %Calculate c*
% c_star = (ao_over_ko) * (rho_o / rho_t) * (a_o / at);