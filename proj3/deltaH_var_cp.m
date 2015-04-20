function [ deltaH ] = deltaH_var_cp( T1, T2, type )

dT=.01;
T=T1;
deltaH=0;

%O2
if (strcmp(type, 'O2'))
    while T < T2
        T = T+dT;
        increment = sp_heats_O2(T) .* dT;
        deltaH = deltaH + increment;
    end
end

%N2
if (strcmp(type, 'N2'))
    while T < T2
        T = T+dT;
        increment = sp_heats_N2(T) .* dT;
        deltaH = deltaH + increment;
    end
end

%CO2
if (strcmp(type, 'CO2'))
    while T < T2
        T = T+dT;
        increment = sp_heats_CO2(T) .* dT;
        deltaH = deltaH + increment;
    end
end

%H2
if (strcmp(type, 'H2O'))
    while T < T2
        T = T+dT;
        increment = sp_heats_H2O(T) .* dT;
        deltaH = deltaH + increment;
    end
end

end