function [ deltaH ] = deltaH_var_cp( T1, T2, n, type, phi, MM)

dT=.01;
T=T1;
deltaH=zeros(1,n);

for i=1:n
    while T(i) < T2(i)
        T(i) = T(i) + dT;
        if strcmp(type, 'JetA')
            increment = sp_heats_JetA(T(i), phi(i), MM) .* dT;
        else
            increment = sp_heats(T(i), type) .* dT;
        end
        deltaH(i) = deltaH(i) + increment;
    end
end
end