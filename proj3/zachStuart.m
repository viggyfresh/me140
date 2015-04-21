function [Ma, To, T, Po_over_P] = zachStuart(Tm, Po, m_dot, A, RF, type)

To = Tm;
R = 286.9;

if strcmp(type, 'air')
    a = 28.11;
    b = 0.1967*10^-2;
    c = 0.4802*10^-5;
    d = -1.966*10^-9;
else
    a = 1.0038;
    b = 0;
    c = 0;
    d = 0;
end

while (1)
    T = Tm - 0.01;
    mfp = (m_dot ./ A) .* sqrt(R .* To) ./ Po;
    while (1)
        I1 = R * (a .* (To - T) + (b / 2) .* (To.^2 - T.^2) + ...
            (c / 3) .* (To.^3 - T.^3) + (d / 4) .* (To.^4 - T.^4));
        I2 = (a .* log (To ./ T) + b .* (To - T) + (c / 2) ...
            .* (To.^2 - T.^2) + (d / 3) .* (To.^3 - T.^3));
        c_p_ave = I1 ./ (To - T);
        Po_over_P = exp(I2);
        lhs = mfp;
        rhs = Inf;
        Ma = -0.001;
        [~, ~, k, ~] = sp_heats(T , type);
        while (abs(lhs - rhs) / lhs > 0.05)
            Ma = Ma + 0.001;
            rhs = Ma .* sqrt(k) .* (1 + (k * R * Ma^2 ...
                / (2 * c_p_ave)))^(0.5) ./ Po_over_P;
        end
        T_guess = Tm ./ (1 + (RF * k * R * Ma^2 / (2 * c_p_ave)));
        if (abs(T - T_guess) / T < 0.01)
            T = T_guess;
            break;
        end
        T = T_guess;
    end
    To_guess = T * (1 + (k * R * Ma^2 / (2 * c_p_ave)));
    if (abs(To - To_guess) / To < 0.01)
        To = To_guess;
        break;
    end
    To = To_guess;
end
end

