function [Ma] = getMach(Tm, Po, m_dot, A, RF)

R = 286.9;
To = Tm;
To_guess = Inf;

a = 28.11;
b = 0.1967*10^-2;
c = 0.4802*10^-5;
d = -1.966*10^-9;

while (1)
    mfp = (m_dot ./ A) .* sqrt(R .* To) ./ Po;
    
    T = Tm;
    T_guess = Inf;
    while (1)
        I1 = R * (a .* (To - T) + (b / 2) .* (To.^2 - T.^2) + ...
            (c / 3) .* (To.^3 - T.^3) + (d / 4) .* (To.^4 - T.^4));
        I2 = (a .* ln (To ./ T) + b .* (To - T) + (c / 2) .* (To.^2 - T.^2) ...
             + (d / 3) .* (To.^3 - T.^3));
        c_p_ave = I1 ./ (To - T);
        Po_over_P = exp(I2);
        lhs = mfp;
        rhs = Inf;
        Ma = -0.01;
        [~, ~, k, ~] = sp_heats(T);
        while (abs(lhs - rhs) / lhs > 0.01)
            Ma = Ma + 0.01;
            rhs = Ma .* sqrt(k) .* (1 + (k * R * Ma^2 / (2 * c_p_ave)))^(0.5) ...
                ./ Po_over_P;
        end
        T_guess = Tm ./ (1 + (RF * k * R * Ma^2 / (2 * c_p_ave)));
        if (abs(T - T_guess) ./ T < 0.01) 
            break;
        end
        T = T_guess;
    end
    
    To_guess = T * (1 + (k * R * Ma^2 / (2 * c_p_ave)));
    if (abs(To - To_guess) ./ To < 0.01)
        break;
    end
    To = To_guess;

end



end

