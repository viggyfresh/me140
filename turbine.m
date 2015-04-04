function [ vars ] = turbine( vars )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 

W_comp = vars.m_dot_core.*vars.c_p.*(vars.T_03-vars.T_013);
W_fan = vars.m_dot.*vars.c_p.*(vars.T_013-vars.T_02);
W_comp_and_fan = W_comp + W_fan;

vars.T_05 = vars.T_04-W_comp_and_fan./(vars.m_dot_core.*vars.c_p);
vars.T_05s = (vars.c_p.*vars.T_04-(vars.c_p.*vars.T_04-vars.c_p.*vars.T_05)./vars.eta_turb)./vars.c_p;


vars.P_05 = vars.P_04.*(vars.T_05s./vars.T_04).^(vars.k./(vars.k-1));

end