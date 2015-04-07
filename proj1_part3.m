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

%Compressor
vars = compressor_var(vars);

%Combustor
vars.P_04 = vars.P_04_over_03.*vars.P_03;
vars.q_dot = deltaH_var_cp(vars.T_03, vars.T_04);

%mdots
vars.m_dot_bp = vars.m_dot.*10./11;
vars.m_dot_core = vars.m_dot./11;

%Turbine
vars = turbine_var(vars);

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

%intial velocity
[~,~,vars.k,~] = sp_heats(vars.T_0_static);
vars.U_0 = vars.Ma.*sqrt(vars.k.*vars.R.*vars.T_0_static);

%Thrust and specific thrust
vars.F_thrust = vars.m_dot_core.*vars.U_8+vars.m_dot_bp.*vars.U_18-vars.m_dot.*vars.U_0;
vars.spec_thrust=vars.F_thrust./vars.m_dot;

%Thrust-specific fuel consumption
vars.m_dot_fuel=vars.q_dot./vars.lhv;
vars.tsfc=vars.m_dot_fuel./vars.F_thrust;

vars




%%%%%%%%%%%%%%Plotting (for part 4)%%%%%%%%%%%%%%%%%%%%

%four plots:
    %core flow for Cruise and SLS
    %bypass flow for Cruise and SLS

    
% see deltaS.m function





