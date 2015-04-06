function [T] = var_cp_comp(T1,T2s,eta_compfan)

dT=.01;
T=T1;
numerator=0;

while T<T2s
    T=T+dT;
    increment = sp_heats(T) .* dT;
    numerator = numerator + increment;
end

target=numerator./eta_compfan;
T=T1;
left=0;

while left<target
    T=T+dT;
    increment=sp_heats(T).*dT;
    left=left+increment;
end
end

