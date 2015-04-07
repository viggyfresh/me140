function [ vars ] = turbine_var( vars )
%This is currently the same as turbine.m, will use as template to for
%variable prameter
% 

W_comp = vars.m_dot_core.*deltaH_var_cp(vars.T_013, vars.T_03);
W_fan = vars.m_dot.*deltaH_var_cp(vars.T_02, vars.T_013);
W_comp_and_fan = W_comp + W_fan;
deltah = -W_comp_and_fan./vars.m_dot_core; %negative so T4 on lower bound
%FIX THIS, FUNCTION HAS PROBLEM, REF T_given_deltaH function
vars.T_05 = T_given_deltaH(vars.T_04, deltah);

vars.k = 1.4; %temporary

vars.T_05s = (vars.c_p.*vars.T_04-(vars.c_p.*vars.T_04-vars.c_p.*vars.T_05)./vars.eta_turb)./vars.c_p;
vars.P_05 = vars.P_04.*(vars.T_05s./vars.T_04).^(vars.k./(vars.k-1));

end

