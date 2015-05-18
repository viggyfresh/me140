function [ ] = graphthethings()

% Cstar and mix ratio for lab data
load('bradycheated.mat');
i1 = start_index;
i2 = final_index;

Po = mean(chamP(i1:i2)) * 1000 + 101325;
D = 0.605; %inches
D = D / 39.370; %meters
At = pi * D^2 / 4;
t = time(i2) - time(i1); %secs
m_dot_fuel = mfuel / 10^3 / t; %kg
m_dot_O2 = mean(m_dot_O2(i1:i2)); 
mdot = m_dot_fuel + m_dot_O2;

cstar_lab = Po / (mdot / At);
mixRatio_lab = m_dot_O2 / m_dot_fuel;

load('final.mat');

figure;
plot(mixRatio, X_frozen_t(iC2H4,:),'g')
hold on;
plot(mixRatio, X_frozen_t(iO2,:),'b')
plot(mixRatio, X_frozen_t(iCO2,:),'r')
plot(mixRatio, X_frozen_t(iH2O,:),'c')
plot(mixRatio, X_frozen_t(iCO,:),'k')
plot(mixRatio, X_frozen_t(iC,:),'--g')
plot(mixRatio, X_frozen_t(iH2,:),'--b')
plot(mixRatio, X_frozen_t(iH,:),'--r')
plot(mixRatio, X_frozen_t(iO,:),'--c')
plot(mixRatio, X_frozen_t(iOH,:),'--k')
xlabel('Mixture Ratio');
ylabel('Mole Fractions at Nozzle Throat');
title('Mixture Ratio vs. Mole Fractions at Nozzle Throat, Frozen');
legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
set(gcf, 'color', 'white');
plotfixer;

% Plot of mole ratios for dissociated case
figure;
plot(mixRatio, X_dissoc_t(iC2H4,:),'g')
hold on;
plot(mixRatio, X_dissoc_t(iO2,:),'b')
plot(mixRatio, X_dissoc_t(iCO2,:),'r')
plot(mixRatio, X_dissoc_t(iH2O,:),'c')
plot(mixRatio, X_dissoc_t(iCO,:),'k')
plot(mixRatio, X_dissoc_t(iC,:),'--g')
plot(mixRatio, X_dissoc_t(iH2,:),'--b')
plot(mixRatio, X_dissoc_t(iH,:),'--r')
plot(mixRatio, X_dissoc_t(iO,:),'--c')
plot(mixRatio, X_dissoc_t(iOH,:),'--k')
xlabel('Mixture Ratio');
ylabel('Mole Fractions at Nozzle Throat');
title('Mixture Ratio vs. Mole Fractions at Nozzle Throat, Chemical Equilbrium');
legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
set(gcf, 'color', 'white');
plotfixer;

% Plot of mole ratios for frozen case
figure;
plot(mixRatio, X_frozen_e(iC2H4,:),'g')
hold on;
plot(mixRatio, X_frozen_e(iO2,:),'b')
plot(mixRatio, X_frozen_e(iCO2,:),'r')
plot(mixRatio, X_frozen_e(iH2O,:),'c')
plot(mixRatio, X_frozen_e(iCO,:),'k')
plot(mixRatio, X_frozen_e(iC,:),'--g')
plot(mixRatio, X_frozen_e(iH2,:),'--b')
plot(mixRatio, X_frozen_e(iH,:),'--r')
plot(mixRatio, X_frozen_e(iO,:),'--c')
plot(mixRatio, X_frozen_e(iOH,:),'--k')
xlabel('Mixture Ratio');
ylabel('Mole Fractions at Nozzle Exit');
title('Mixture Ratio vs. Mole Fractions at Nozzle Exit, Frozen');
legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
set(gcf, 'color', 'white');
plotfixer;

% Plot of mole ratios for dissociated case
figure;
plot(mixRatio, X_dissoc_e(iC2H4,:),'g')
hold on;
plot(mixRatio, X_dissoc_e(iO2,:),'b')
plot(mixRatio, X_dissoc_e(iCO2,:),'r')
plot(mixRatio, X_dissoc_e(iH2O,:),'c')
plot(mixRatio, X_dissoc_e(iCO,:),'k')
plot(mixRatio, X_dissoc_e(iC,:),'--g')
plot(mixRatio, X_dissoc_e(iH2,:),'--b')
plot(mixRatio, X_dissoc_e(iH,:),'--r')
plot(mixRatio, X_dissoc_e(iO,:),'--c')
plot(mixRatio, X_dissoc_e(iOH,:),'--k')
xlabel('Mixture Ratio');
ylabel('Mole Fractions at Nozzle Exit');
title('Mixture Ratio vs. Mole Fractions at Nozzle Exit, Chemical Equilbrium');
legend('C_2H_4', 'O_2', 'CO_2', 'H_2O', 'CO', 'C', 'H_2', 'H', 'O', 'OH');
set(gcf, 'color', 'white');
plotfixer;

% Plot of Throat temperature and stag temperature
figure;
plot(mixRatio, To, 'r', mixRatio, T_t_frozen, '--b', mixRatio, T_t_dissoc, 'b', mixRatio, T_e_frozen, '--g', mixRatio, T_e_dissoc, 'g');
xlabel('Mixture Ratio');
ylabel('Temperature (K)');
title('Mixture Ratio vs. Various Temperatures');
legend('T_0', 'T_t frozen', 'T_t', 'T_e frozen', 'T_e');
set(gcf, 'color', 'white');
plotfixer;

% Plot of c star
figure;
plot(mixRatio, c_star_frozen, '--k', mixRatio, c_star_dissoc, 'k');
hold on;
plot(mixRatio_lab, cstar_lab, '*', 'markersize', 25);
xlabel('Mixture Ratio');
ylabel('c^* (m/s)');
title('c^* vs. Mixture Ratio');
legend('c^* Frozen', 'c^*', 'c^* stock motor run');
ylim([0 6000])
set(gcf, 'color', 'white');
plotfixer;
yL = get(gca, 'YLim');
line([mixRatio_lab, mixRatio_lab], yL, 'Linestyle', '--');
plotfixer;

hold off;

%Plot of Velocity
figure;
plot(mixRatio, V_e_frozen, '--m', mixRatio, V_e_dissoc, 'm');
title('Exit Velocity vs. Mixture Ratio');
legend('V_e frozen', 'V_e');
set(gcf, 'color', 'white');
plotfixer;

%Plot thrust coefficient
figure;
plot(mixRatio, Cf_dissoc);
xlabel('Mixture Ratio');
ylabel('Thrust Coefficient');
title('Thrust Coefficient vs. Mixture Ratio');
set(gcf, 'color', 'white');
plotfixer;

%Plot optimal nozzle expansion ratio
figure;
plot(mixRatio, epsilon_dissoc);
xlabel('Mixture Ratio');
ylabel('Ratio');
title('Optimal Nozzle Expansion Ratio');
set(gcf, 'color', 'white');
plotfixer;

%Plot everything
%Plot of Throat temperature and stag temperature
figure;
plot(mixRatio, To, 'r', mixRatio, T_t_frozen, '--b', mixRatio, T_t_dissoc, 'b', mixRatio, T_e_frozen, '--g', mixRatio, T_e_dissoc, 'g');
xlabel('Mixture Ratio');
ylabel('Temperature (K), Speed (m/s)');
title('Mixture Ratio vs. Various Quantities');
%legend('T_0', 'T_t frozen', 'T_t', 'T_e frozen', 'T_e);
set(gcf, 'color', 'white');
% Plot of c star
hold on;
plot(mixRatio, c_star_frozen, '--k', mixRatio, c_star_dissoc, 'k');
% legend('c^* Frozen', 'c^*');
ylim([0 6000])
set(gcf, 'color', 'white');
%Plot of Velocity
hold on;
plot(mixRatio, V_e_frozen, '--m', mixRatio, V_e_dissoc, 'm');

legend('T_0', 'T_t frozen', 'T_t', 'T_e frozen', 'T_e', ...
    'C^* frozen', 'c^*', 'V_e frozen', 'V_e');

plotfixer;

end

