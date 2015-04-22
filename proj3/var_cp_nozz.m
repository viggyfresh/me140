function [T] = var_cp_nozz(T1,T2s,eta_nozz)

dT=.01;
T=T2s;
denominator=[0 0];

for i=1:1
    while T(i)<T1(i)
        T(i)=T(i)+dT;
        increment = sp_heats(T(i), 'air') .* dT;
        denominator(i) = denominator(i) + increment;
    end

    target(i)=-denominator(i).*eta_nozz(i);
    T(i)=T1(i);
    left=0;

    while left>target(i)
        T(i)=T(i)-dT;
        increment=sp_heats(T(i), 'air')*dT;
        left=left-increment;
    end
end
end

