function [ ] = graphthethings( file )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load(file);

% Plot of mole ratios for frozen case
% figure;
% plot(mixRatio, X_frozen_t(iC2H4,:),'g')
% hold on;
% plot(mixRatio, X_frozen_t(iO2,:),'b')
% plot(mixRatio, X_frozen_t(iCO2,:),'r')
% plot(mixRatio, X_frozen_t(iH2O,:),'c')
% plot(mixRatio, X_frozen_t(iCO,:),'k')
% plot(mixRatio, X_frozen_t(iC,:),'--g')
% plot(mixRatio, X_frozen_t(iH2,:),'--b')
% plot(mixRatio, X_frozen_t(iH,:),'--r')
% plot(mixRatio, X_frozen_t(iO,:),'--c')
% plot(mixRatio, X_frozen_t(iOH,:),'--k')
% xlabel('Mixture Ratio');
% ylabel('Mole Fraction at Nozzle Throat');
% title('Mixture Ratio vs. Mole Fractions at Nozzle Throat, Frozen');
% legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
% set(gcf, 'color', 'white');git ad
% plotfixer;
% 
% % Plot of mole ratios for dissociated case
% figure;
% plot(mixRatio, X_dissoc_t(iC2H4,:),'g')
% hold on;
% plot(mixRatio, X_dissoc_t(iO2,:),'b')
% plot(mixRatio, X_dissoc_t(iCO2,:),'r')
% plot(mixRatio, X_dissoc_t(iH2O,:),'c')
% plot(mixRatio, X_dissoc_t(iCO,:),'k')
% plot(mixRatio, X_dissoc_t(iC,:),'--g')
% plot(mixRatio, X_dissoc_t(iH2,:),'--b')
% plot(mixRatio, X_dissoc_t(iH,:),'--r')
% plot(mixRatio, X_dissoc_t(iO,:),'--c')
% plot(mixRatio, X_dissoc_t(iOH,:),'--k')
% xlabel('Mixture Ratio');
% ylabel('Mole Fraction at Nozzle Throat');
% title('Mixture Ratio vs. Mole Fractions at Nozzle Throat, Chemical Equilbrium');
% legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
% set(gcf, 'color', 'white');
% plotfixer;
% 
% % Plot of mole ratios for frozen case
% figure;
% plot(mixRatio, X_frozen_e(iC2H4,:),'g')
% hold on;
% plot(mixRatio, X_frozen_e(iO2,:),'b')
% plot(mixRatio, X_frozen_e(iCO2,:),'r')
% plot(mixRatio, X_frozen_e(iH2O,:),'c')
% plot(mixRatio, X_frozen_e(iCO,:),'k')
% plot(mixRatio, X_frozen_e(iC,:),'--g')
% plot(mixRatio, X_frozen_e(iH2,:),'--b')
% plot(mixRatio, X_frozen_e(iH,:),'--r')
% plot(mixRatio, X_frozen_e(iO,:),'--c')
% plot(mixRatio, X_frozen_e(iOH,:),'--k')
% xlabel('Mixture Ratio');
% ylabel('Mole Fraction at Nozzle Exit');
% title('Mixture Ratio vs. Mole Fractions at Nozzle Exit, Frozen');
% legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
% set(gcf, 'color', 'white');
% plotfixer;
% 
% % Plot of mole ratios for dissociated case
% figure;
% plot(mixRatio, X_dissoc_e(iC2H4,:),'g')
% hold on;
% plot(mixRatio, X_dissoc_e(iO2,:),'b')
% plot(mixRatio, X_dissoc_e(iCO2,:),'r')
% plot(mixRatio, X_dissoc_e(iH2O,:),'c')
% plot(mixRatio, X_dissoc_e(iCO,:),'k')
% plot(mixRatio, X_dissoc_e(iC,:),'--g')
% plot(mixRatio, X_dissoc_e(iH2,:),'--b')
% plot(mixRatio, X_dissoc_e(iH,:),'--r')
% plot(mixRatio, X_dissoc_e(iO,:),'--c')
% plot(mixRatio, X_dissoc_e(iOH,:),'--k')
% xlabel('Mixture Ratio');
% ylabel('Mole Fraction at Nozzle Exit');
% title('Mixture Ratio vs. Mole Fractions at Nozzle Exit, Chemical Equilbrium');
% legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
% set(gcf, 'color', 'white');
% plotfixer;

% Plot of Throat temperature and stag temperature
figure;
plot(mixRatio, To, 'r', mixRatio, T_t_frozen, '--b', mixRatio, T_t_dissoc, 'b', mixRatio, T_e_frozen, '--g', mixRatio, T_e_dissoc, 'g');
xlabel('Mixture Ratio');
ylabel('Temperature (K)');
title('Mixture Ratio vs. Various Temperatures');
%legend('T_0', 'T_t frozen', 'T_t', 'T_e frozen', 'T_e);
set(gcf, 'color', 'white');
plotfixer;

% Plot of c star
hold on;
plot(mixRatio, c_star_frozen, '--k', mixRatio, c_star_dissoc, 'k');
xlabel('Mixture Ratio');
ylabel('c^* (m/s)');
title('c^* vs. Mixture Ratio');
% legend('c^* Frozen', 'c^*');
ylim([0 6000])
set(gcf, 'color', 'white');

%Plot of Velocity
hold on;
plot(mixRatio, V_e_frozen, '--m', mixRatio, V_e_dissoc, 'm');

legend('T_0', 'T_t frozen', 'T_t', 'T_e frozen', 'T_e', ...
    'C^* frozen', 'c^*', 'V_e frozen', 'V_e');

plotfixer;

%Plot thrust coefficient
figure;
plot(mixRatio, Cf_dissoc);
xlabel('Mixture Ratio');
ylabel('Thrust Coefficient');
title('Thrust Coefficient vs. Mixture Ratio');
plotfixer;

%Plot optimal nozzle expansion ratio
figure;
plot(mixRatio, epsilon_dissoc);
xlabel('Mixture Ratio');
ylabel('Ratio');
title('Optimal Nozzle Expansion Ratio');
plotfixer;


end

