function [cp] = cpGuess(T1,P1,P2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

R=287;
cp_guess=1000;
[cp1,~,~,~]=sp_heats(T1);
cp1=cp1*1000
cp_av=0.5*(cp_guess+cp1)


cp=0;
end

