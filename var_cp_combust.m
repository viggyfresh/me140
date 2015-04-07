function [ deltaH ] = var_cp_combust( T1,T2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dT=.01;
T=T1;
deltaH=0;

while T<T2
    T=T+dT;
    increment = sp_heats(T) .* dT;
    deltaH = deltaH + increment;
end

end

