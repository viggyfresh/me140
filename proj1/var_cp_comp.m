function [T] = var_cp_comp(T1,T2s,eta_compfan)

dT=.01;
T=T1;
numerator=[0 0];

for i=1:2
    while T(i)<T2s(i)
        T(i)=T(i)+dT;
        increment = sp_heats(T(i)) .* dT;
        numerator(i) = numerator(i) + increment;
    end

    target(i)=numerator(i)./eta_compfan(i);
    T(i)=T1(i);
    left=0;

    while left<target(i)
        T(i)=T(i)+dT;
        increment=sp_heats(T(i)).*dT;
        left=left+increment;
    end
end
end

