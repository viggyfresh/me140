function [ T ] = T_given_deltaH(T1, deltah)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dT = 0.01;
T=T1;
left=0;
target = deltah;

%PROBLEM: Turbine starts at higher temp than exit = negative number
while abs(left-target)>dT
    T=T+dT;
    increment=sp_heats(T).*dT;
    left=left+increment;
end

end

