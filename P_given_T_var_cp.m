function [P2] = P_given_T_var_cp(T1, T2, P1)
%T2 < T1 for turbine

R = 286.9;
T = T2;
dT = 0.01;
cp_int = [0 0];


for i=1:2
    while T(i) < T1(i)
        T(i) = T(i) + dT;
        increment = sp_heats(T(i))*(dT/T(i));
        cp_int(i) = cp_int(i) + increment;
    end
    P2(i) = P1(i)/(exp(cp_int(i)/R));
end

end

