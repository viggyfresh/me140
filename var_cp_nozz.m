function [T] = var_cp_nozz(T1,T2s,eta_nozz)

dT=.01;
T=T2s; %because T2<T1
denominator=[0 0];

for i=1:2
    while T(i)<T1(i)
        T(i)=T(i)+dT;
        increment = sp_heats(T(i)) .* dT;
        denominator(i) = denominator(i) + increment;
    end

    target(i)=-denominator(i).*eta_nozz(i);
    T(i)=T1(i);
    left=0;

    % We are looking for lower bound, which is why target is negative
    while left>target(i)
        T(i)=T(i)-dT;
        increment=sp_heats(T(i))*dT;
        left=left-increment;
    end
end
end

