function [To2s] = comp_Ts(To1,Po2_over_Po1, n)

R = 286.9;
target=R.*log(Po2_over_Po1);

T2=To1;
dT=0.01;
cp_int=zeros([1, n]);

for i=1:n
    while (cp_int(i) < target(i))
        T2(i) = T2(i) + dT;
        dcp_int = sp_heats(T2(i))*(dT/T2(i));
        cp_int(i) = cp_int(i) + dcp_int;
    end 
end
To2s = T2;
end

