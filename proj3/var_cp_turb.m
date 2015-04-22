function [T] = var_cp_turb(T1,T2,eta_turb)

dT=.01;
T=T2;
numerator=[0 0];

for i=1:1
    while T(i)<T1(i)
        T(i)=T(i)+dT;
        increment = sp_heats(T(i), 'air') .* dT;
        numerator(i) = numerator(i) + increment;
    end

    target(i)=-eta_turb(i)/numerator(i);
    T(i)=T1(i);
    left=0;

    while left>target(i)
        T(i)=T(i)-dT;
        increment=sp_heats(T(i), 'air')*dT;
        left=left-increment;
    end
end
end

