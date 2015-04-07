function [ T ] = T_given_deltaH(T1, deltah)
%T2 < T1
dT = 0.01;
T=T1;
left=[0 0];
target = deltah;

%PROBLEM: Turbine starts at higher temp than exit = negative number
for i = 1:2
    while left(i) > target(i)
        T(i)=T(i)-dT;
        increment=sp_heats(T(i)).*dT;
        left(i)=left(i)-increment;
        print = left(i) - target(i);
    end
end

end
