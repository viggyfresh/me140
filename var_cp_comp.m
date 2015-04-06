function [T] = var_cp(T1,T2s,eta_compfan)

R = 286.9;
cp_2_calc = 0;
cp_2_guess = 500;
cp_avg_s = 1/2 .* ( sp_heats(T1) + sp_heats(T2s) );

count = 0;
while abs(cp_2_calc-cp_2_guess)>0.01
    cp_2_guess = cp_2_calc;
    cp_avg = 1/2 .* ( sp_heats(T1) + cp_2_guess );
    T2_calc = T1 + (cp_avg_s.*(T2s-T1))./(cp_avg.*eta_compfan);
    cp_2_calc = sp_heats(T2_calc);
    count = count+1;
end
cp = cp_2_calc;
T = T2_calc;
end

