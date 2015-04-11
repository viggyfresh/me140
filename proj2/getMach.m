function [Ma] = getMach(Tm, Po, m_dot, A, RF)


%%OLD CODE

To = Tm;
To_guess = Inf;

R = 286.9;
a = 28.11;
b = 0.1967*10^-2;
c = 0.4802*10^-5;
d = -1.966*10^-9;

while (1)    
    T = 1;
    T_guess = Inf;
    while (1)
        mfp = (m_dot ./ A) .* sqrt(R .* To) ./ Po;
        I1 = R * (a .* (To - T) + (b / 2) .* (To.^2 - T.^2) + ...
            (c / 3) .* (To.^3 - T.^3) + (d / 4) .* (To.^4 - T.^4));
        I2 = (a .* log (To ./ T) + b .* (To - T) + (c / 2) .* (To.^2 - T.^2) ...
             + (d / 3) .* (To.^3 - T.^3));
        % what do we use instead of c_p_ave?
        c_p_ave = I1 ./ (To - T);
        Po_over_P = exp(I2);
        lhs = mfp
        rhs = Inf;
        Ma = -0.01;
        [~, ~, k, ~] = sp_heats(T);
        while (abs(lhs - rhs) / lhs > 0.01 && Ma < 1)
            Ma = Ma + 0.01;
        % rhs = MFP
            [Po_over_P2, To_over_T, rhs] = the_var(Ma, T);
            rhs;
%             rhs = Ma .* sqrt(k) .* (1 + (k * R * Ma^2 / (2 * c_p_ave)))^(0.5) ...
%                   ./ Po_over_P;
             
            
        end
        rhs
%         if (Ma > 1)
%             T = T + 100;
%             continue;
%         end
        T_guess = Tm ./ (1 + (RF * k * R * Ma^2 / (2 * c_p_ave))); %%might need to remove cp_avg
        %To_over_T = (Tm ./ T -1)./RF+1; %don't know if this works
        %T_guess = To / To_over_T;
        if (abs(T - T_guess) / T < 0.01) 
            break;
        end
        T = T_guess;
    end
    
    To_guess = T * (1 + (k * R * Ma^2 / (2 * c_p_ave))); %%might need to remove cp_avg
    %To_guess = T .* To_over_T;
    if (abs(To - To_guess) / To < 0.01)
        break;
    end
    To = To_guess;

end
end

