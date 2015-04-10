function [P_stag_ratio, T_stag_ratio, MFP] = the_var(Ma, T)
[cp1, cv1, gamma1] = sp_heats(T);

R = 287;
c = sqrt(gamma1*R*T);
U=Ma*c;

target = 0.5*U^2;
cp_int = 0; %integral needed to find cp
pressure_int = 0; %integral needed to find pressure ratio
To = T; %To must be at least T
dT = 0.01; %Increment size

% Increment until To is found
while (cp_int < target)
    To = To + dT;
    dcp_int = sp_heats(To)*dT;
    cp_int = cp_int + dcp_int;
    pressure_int = pressure_int + dcp_int/(R*To);
end 
P_stag_ratio= exp(pressure_int);
T_stag_ratio = To/T;
MFP = Ma*sqrt(gamma1)*sqrt(T_stag_ratio)/P_stag_ratio;
end

