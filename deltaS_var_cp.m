function [deltaS] = deltaS_var_cp(T1, T2, P1, P2)

R = 286.9;
T = T1;
dT = 0.01;
cp_int = [0 0];

%Do we need to edit this function for the case that T2 is smaller than T1?
%The while loop below wont trip for that case as is (I don't think)

for i=1:2
    
    if T1(i)<T2(i)
        while T(i) < T2(i)
            T(i) = T(i) + dT;
            increment = sp_heats(T(i))*(dT/T(i));
            cp_int(i) = cp_int(i) + increment;
        end
        deltaS(i) = cp_int(i) - R*log(P2(i)/P1(i));
    else
        while T(i) > T2(i)
            T(i) = T(i) - dT;
            increment = sp_heats(T(i))*(dT/T(i));
            cp_int(i) = cp_int(i) - increment;
        end
        deltaS(i) = cp_int(i) - R*log(P2(i)/P1(i));
    end
end

end

