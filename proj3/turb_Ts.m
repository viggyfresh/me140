function [To2s] = turb_Ts(To1,Po2_over_Po1, n, type, phi, MM)

R = 286.9;
target=R*log(Po2_over_Po1);

T2=To1;
dT=0.01;
cp_int=zeros([1,n]);
for i=1:n
    while (cp_int(i) > target(i))
        T2(i) = T2(i) - dT;
        if strcmp(type, 'JetA')
            dcp_int = sp_heats_JetA(T2(i), phi(i), MM) .* (dT./T2(i));
        else
            dcp_int = sp_heats(T2(i), type).*(dT./T2(i));
        end
        cp_int(i) = cp_int(i) - dcp_int;
    end 
end

To2s = T2;
end

