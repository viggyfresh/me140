function [ vars ] = compressor_var( vars )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 

P03_over_P02 = vars.P_ratio_overall;
P013_over_P02 = vars.P_ratio_fan;
P03_over_P013 = P03_over_P02 ./ P013_over_P02;
vars.P_03 = vars.P_013 .* P03_over_P013;

vars.T_03s = var_cp(vars.T_013,P03_over_P013);
vars.T_03 = var_cp_comp(vars.T_013,vars.T_03s,vars.eta_comp);

end