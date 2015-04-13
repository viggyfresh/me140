function [ deltaH ] = deltaH_var_cp( T1,T2,n )

dT=.01;
T=T1;
deltaH=zeros(1,n);

for i= 1:n
    while T(i)<T2(i)
        T(i)=T(i)+dT;
        increment = sp_heats(T(i)) .* dT;
        deltaH(i) = deltaH(i) + increment;
    end
end

end