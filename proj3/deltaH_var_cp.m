function [ deltaH ] = deltaH_var_cp( T1, T2, type )

dT=.01;
T=T1;
deltaH=0;

while T < T2
    T = T+dT;
    increment = sp_heats(T,type) .* dT;
    deltaH = deltaH + increment;
end

end