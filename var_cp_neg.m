function [T] = var_cp_neg(T1,P2_over_P1)
%when P2<P1
R = 286.9;
cp_2_calc = 0;
cp_2_guess = 500;
count = 0;
target=R*log(P2_over_P1);

T2=T1;
dT=0.01;
cp_int=[0 0];
for i=1:2
    while (cp_int(i) > target(i))
        T2(i) = T2(i) - dT;
        dcp_int = sp_heats(T2(i)).*(dT./T2(i));
        cp_int(i) = cp_int(i) - dcp_int;
    end 
end

T = T2;
end

