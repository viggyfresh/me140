function [Ma] = getMach(Tm,pressure,airFlow,area)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

R=286.9;
To=Tm;
mfp=(airFlow./area).*sqrt((R*To)./pressure);


end

