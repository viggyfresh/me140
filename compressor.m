function [ vars ] = compressor( vars )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 

P03_over_P02 = vars.P_ratio_overall;
P013_over_P02 = vars.P_ratio_fan;
P03_over_P013 = P03_over_P02 ./ P013_over_P02;
vars.P_03 = vars.P_013 .* P03_over_P013;

T_03s = vars.T_013 .* P03_over_P013 .^ ((vars.k - 1)./vars.k);
vars.T_03 = ((vars.c_p .* T_03s - vars.c_p .* vars.T_013) ./ vars.eta_comp + (vars.c_p * vars.T_013)) ./ vars.c_p;

end