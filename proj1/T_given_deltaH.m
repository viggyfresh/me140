function [ T ] = T_given_deltaH(T1, deltah)
dT = 0.01;
T=T1;
left=[0 0];
target = deltah;

for i = 1:2
    while left(i) > target(i)
        T(i)=T(i)-dT;
        increment=sp_heats(T(i)).*dT;
        left(i)=left(i)-increment;
        print = left(i) - target(i);
    end
end

end
