function [Ma, To, T, Po_over_P] = viggyFresh(Tm, Po, m_dot, A, RF, phi, MM)

[~, ~, ~, R] = sp_heats_JetA(Tm, phi, MM);
To = Tm;
dT = 0.01;

while (1)
    T = Tm - 0.01;
    mfp = (m_dot ./ A) .* sqrt(R .* To) ./ Po;
    while (1)
        I1 = 0;
        I2 = 0;
        T_prime = T;
        dT = 0.01;
        while T_prime < To
            T_prime = T_prime + dT;
            [cp, ~, ~, R] = sp_heats_JetA(T_prime, phi, MM);
            I1 = I1 + (cp * dT);
            I2 = I2 + (cp / (R * T_prime)) * dT;
        end
        c_p_ave = I1 ./ (To - T);
        Po_over_P = exp(I2);
        lhs = mfp;
        rhs = Inf;
        Ma = -0.001;
        [~, ~, k, ~] = sp_heats_JetA(T, phi, MM);
        while (abs(lhs - rhs) / lhs > 0.07)
            Ma = Ma + 0.001;
            %[~, ~, rhs] = the_var_JetA(Ma, T, phi, MM);
            rhs = Ma .* sqrt(k) .* (1 + (k * R * Ma^2 ...
                  / (2 * c_p_ave)))^(0.5) ./ Po_over_P;
            abs(lhs-rhs) / lhs
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



