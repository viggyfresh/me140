clc;
clear all;
close all;

% Want: thrust, specific thrust, thurst specific fuel consumption

% Table 1: Parameters for cruise (index 1) and SLS (index 2)

vars.h = [10.67 * 10^3 0]; % m
vars.T_0_static = [218.8 288.15]; % K
vars.P_0_static = [0.239 1.014] * 10^5; % Pa
vars.Ma = [0.78 0];
vars.P_ratio_overall = [32 28];
vars.P_ratio_fan = [1.55 1.52];
vars.T_04 = [1450 1650]; % K
vars.m_dot = [110 265]; % kg/s
vars.BR = [10 10];

% Table 2: Engine Data

vars.P_02_over_00 = [0.998 1.00];
vars.eta_fan = [0.95 0.95];
vars.eta_comp = [0.89 0.89];
vars.eta_turb = [0.90 0.90];
vars.eta_nozz = [0.95 0.95];
vars.P_04_over_03 = [0.95 0.95];

%Other useful stuff
vars.c_p = 1005; % J / kg * K
vars.R = 287;
vars.lhv= 42.8 * 10^6; %J/kg

%Pre fan and compressor calculations
for i=1:2
    [pressRatio_0_static(i),tempRatio_0_static(i),~]=the_var(vars.Ma(i),...
        vars.T_0_static(i));
end
vars.P_00=vars.P_0_static.*pressRatio_0_static;
vars.P_02=vars.P_02_over_00.*vars.P_00;
vars.T_00=vars.T_0_static.*tempRatio_0_static;
vars.T_02=vars.T_00;

[vars.T_013s]=var_cp(vars.T_02,vars.P_ratio_fan);
vars.T_013=var_cp_comp(vars.T_02,vars.T_013s,vars.eta_fan);
vars.P_013=vars.P_02.*vars.P_ratio_fan;

%Calculate entropy change from 0 to 2
vars.deltaS_0to2=deltaS_var_cp(vars.T_00,vars.T_02,vars.P_00,vars.P_02);

%Calculate entropy change from 2 to 13
vars.deltaS_2to13=deltaS_var_cp(vars.T_02,vars.T_013,vars.P_02,vars.P_013);

%Compressor
vars = compressor_var(vars);

%Calculate entropy change from 13 to 3
vars.deltaS_13to3=deltaS_var_cp(vars.T_013,vars.T_03,vars.P_013,vars.P_03);

%mdots
vars.m_dot_bp = vars.m_dot.*10./11;
vars.m_dot_core = vars.m_dot./11;

%Combustor
vars.P_04 = vars.P_04_over_03.*vars.P_03;
vars.q_dot = vars.m_dot_core.*deltaH_var_cp(vars.T_03, vars.T_04);

%Calculate entropy change from 3 to 4
vars.deltaS_3to4=deltaS_var_cp(vars.T_03,vars.T_04,vars.P_03,vars.P_04);


%Turbine
vars = turbine_var(vars);

%Calculate entropy change from 4 to 5
vars.deltaS_4to5=deltaS_var_cp(vars.T_04,vars.T_05,vars.P_04,vars.P_05);

%Core Nozzle
vars.T_07 = vars.T_05;
vars.P_07 = vars.P_05;
vars.P_8 = vars.P_0_static;
vars.T_8s = var_cp_neg(vars.T_07, vars.P_8./vars.P_07);
vars.T_8 = var_cp_nozz(vars.T_07, vars.T_8s, vars.eta_nozz);
vars.U_8 = sqrt(2.*deltaH_var_cp(vars.T_8, vars.T_07));

%BP Nozzle
vars.P_18 = vars.P_0_static;
vars.T_18s = var_cp_neg(vars.T_013, vars.P_18./vars.P_013);
vars.T_18 = var_cp_nozz(vars.T_013, vars.T_18s, vars.eta_nozz);
vars.U_18 = sqrt(2.*deltaH_var_cp(vars.T_18, vars.T_013));

%Calculate entropy change from 7 to 8
vars.deltaS_5to8=deltaS_var_cp(vars.T_07,vars.T_8,vars.P_07,vars.P_8);

%Calculate entrop change from 13 to 18
vars.deltaS_13to18=deltaS_var_cp(vars.T_013,vars.T_18,vars.P_013,vars.P_18);

%intial velocity
[~,~,vars.k,~] = sp_heats(vars.T_0_static);
vars.U_0 = vars.Ma.*sqrt(vars.k.*vars.R.*vars.T_0_static);

%Thrust and specific thrust
vars.F_thrust = vars.m_dot_core.*vars.U_8+vars.m_dot_bp.*vars.U_18-vars.m_dot.*vars.U_0;
vars.spec_thrust=vars.F_thrust./vars.m_dot;

%Thrust-specific fuel consumption
vars.m_dot_fuel=vars.q_dot./vars.lhv;
vars.tsfc=vars.m_dot_fuel./vars.F_thrust;

%%%%%%%%%%%%%%Plotting (for part 4)%%%%%%%%%%%%%%%%%%%%

entropy_state0=[0 0];
entropy_state2=entropy_state0+vars.deltaS_0to2;
entropy_state13=entropy_state2+vars.deltaS_2to13;
entropy_state3=entropy_state13+vars.deltaS_13to3;
entropy_state4=entropy_state3+vars.deltaS_3to4;
entropy_state5=entropy_state4+vars.deltaS_4to5;
entropy_state8=entropy_state5+vars.deltaS_5to8;


% Keep bypass separate
entropy_state18=entropy_state13+vars.deltaS_13to18;

%Vectors
vars.entropy_states=[entropy_state0; entropy_state2; entropy_state13; ...
    entropy_state3; entropy_state4; entropy_state5; entropy_state8];

vars.temp_states=[vars.T_0_static; vars.T_02; vars.T_013; vars.T_03;...
    vars.T_04; vars.T_05; vars.T_8];

vars.entropy_bp = [entropy_state0; entropy_state2; entropy_state13; entropy_state18];
vars.temp_bp = [vars.T_0_static; vars.T_02; vars.T_013; vars.T_18];

vars

%Plot for cruise core
figure;
labels = {'0'; '02'; '013'; '03'; '04'; '05'; '8'};
plot(vars.entropy_states(:,1), vars.temp_states(:,1),'LineStyle', '--','marker','.','Markersize',20,'color', 'k');
text(vars.entropy_states(:,1), vars.temp_states(:,1), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'Fontsize', 14);
xlabel('deltaS (referenced from state 0) (J/kg*K)', 'FontSize',14);
ylabel('Temperature (K)', 'FontSize',14);
title ('T-S Graph for Cruise Core Flow', 'FontSize',14);
set(gcf, 'color', 'white');

%Plot for cruise BP
figure;
labels = {'0'; '02'; '013'; '18'};
plot(vars.entropy_bp(:,1), vars.temp_bp(:,1),'LineStyle', '--','marker','.','Markersize',20,'color', 'k');
text(vars.entropy_bp(:,1), vars.temp_bp(:,1), labels, 'VerticalAlignment','bottom','HorizontalAlignment','left','Fontsize', 14);
xlabel('deltaS (referenced from state 0) (J/kg*K)','FontSize',14);
ylabel('Temperature (K)','FontSize',14);
title ('T-S Graph for Cruise Bypass Flow','FontSize',14);
set(gcf, 'color', 'white');

%Plot for SLS core
figure;
labels = {'0'; '02'; '013'; '03'; '04'; '05'; '8'};
plot(vars.entropy_states(:,2), vars.temp_states(:,2), 'LineStyle', '--', 'marker','.','Markersize',20,'color', 'r');
text(vars.entropy_states(:,2), vars.temp_states(:,2), labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'Fontsize', 14);
xlabel('deltaS (referenced from state 0) (J/kg*K)','FontSize',14);
ylabel('Temperature (K)','FontSize',14);
title ('T-S Graph for SLS Core Flow','FontSize',14);
set(gcf, 'color', 'white');

%Plot for SLS BP
figure;
labels = {'0'; '02'; '013'; '18'};
plot(vars.entropy_bp(:,2), vars.temp_bp(:,2), 'LineStyle', '--', 'marker','.','Markersize',20,'color', 'r');
text(vars.entropy_bp(:,2), vars.temp_bp(:,2), labels, 'VerticalAlignment','bottom','HorizontalAlignment','left','Fontsize', 14);
xlabel('deltaS (referenced from state 0) (J/kg*K)','FontSize',14);
ylabel('Temperature (K)','FontSize',14);
title ('T-S Graph for SLS Bypass Flow','FontSize',14);
set(gcf, 'color', 'white');





