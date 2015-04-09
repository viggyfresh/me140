function [ vars ] = turbine_var( vars )

W_comp = vars.m_dot_core.*deltaH_var_cp(vars.T_013, vars.T_03);
W_fan = vars.m_dot.*deltaH_var_cp(vars.T_02, vars.T_013);
W_comp_and_fan = W_comp + W_fan;
deltah = -W_comp_and_fan./vars.m_dot_core; 

vars.T_05 = T_given_deltaH(vars.T_04, deltah);


vars.T_05s = var_cp_turb(vars.T_04, vars.T_05, vars.eta_turb);

vars.P_05 = P_given_T_var_cp(vars.T_04, vars.T_05s, vars.P_04);


end