function [T] = var_cp_neg(T1,P2_over_P1)

R = 286.9;
cp_2_calc = 0;
cp_2_guess = 500;
count = 0;
target=R*log(P2_over_P1);

T2=T1;
dT=0.01;
cp_int=0;
while (cp_int > target)
    T2 = T2 - dT;
    dcp_int = sp_heats(T2, 'air').*(dT./T2);
    cp_int = cp_int - dcp_int;
end


T = T2;
end

