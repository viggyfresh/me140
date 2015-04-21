function [Ma, T, P_ratio] = richieTran(To, Po, mdot, A)

R = 286.9; %J/kg*K

mfp = (mdot ./ A) .* sqrt(R .* To) ./ Po; %Target

Ma = 0;
dMa = 0.001;
mfp_guess = Inf;
toler = 0.05;
T_ratio_guess = 1;

while (abs(mfp - mfp_guess) / mfp > toler)
    Ma = Ma + dMa; %increment Mach
    To_guess = 0;
    while (abs (To - To_guess) / To > toler) %iterate throuh temperatures
        T = To / T_ratio_guess;
        [P_ratio, T_ratio_guess, mfp_guess] = the_var(Ma, T);
        To_guess = T_ratio_guess * T;
    end
end
end
