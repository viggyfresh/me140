function [mfp, Ma, T] = richieTran( To, Po, mdot, A, n)

R = 286.9; %J/kg*K

mfp = (mdot ./ A) .* sqrt(R .* To) ./ Po; %Target

for i=1:n
    Ma = zeros(1,n);
    dMa = 0.01;
    mfp_guess = 0;
    toler = 0.01;
    T_ratio_guess = 1; 

    while (abs(mfp(i) - mfp_guess) / mfp(i) > toler)
        Ma(i) = Ma(i) + dMa; %increment Mach
        To_guess = 0;
        while (abs (To(i) - To_guess(i)) / To(i) > toler) %iterate throuh temperatures
            T(i) = To(i) / T_ratio_guess;
            [P_ratio_guess, T_ratio_guess, mfp_guess] = the_var(Ma(i), T(i));
            To_guess = T_ratio_guess * T(i);
        end
    end
end
end

