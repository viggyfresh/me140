function [ deltaH ] = deltaH_var_cp( T1,T2 )

dT=.01;
T=T1;
deltaH=[0 0];

for i=1:2
    while T(i)<T2(i)
        T(i)=T(i)+dT;
        increment = sp_heats(T(i)) .* dT;
        deltaH(i) = deltaH(i) + increment;
    end
end

end