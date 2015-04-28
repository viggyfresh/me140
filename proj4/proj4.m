clc;
clear all;
close all;

% Part 1
alpha = 0;
lambda = 2;
P_standard = 101.325 * 10^3;
P = P_standard;
T_standard = 298;
T_series = 25:5:1000;
T_series = T_series + 273;

for i=1:length(T_series)
    T = T_series(i);
    [eta.LHV(i), eta.HHV(i), eta.actual(i)] = lucio(T, P, alpha, lambda);
    eta.carnot(i) = (1 - (T_standard / T));
end

% Something (might be) wrong with these values
% Plotting for Part 1
figure;
plot(T_series, eta.LHV * 100, T_series, eta.HHV * 100,...
     T_series, eta.actual * 100, T_series, eta.carnot * 100);
xlabel('Temperature (K)');
ylabel('Efficiency (%)');
legend('LHV', 'HHV', '\DeltaH', '\eta_c');
title('First Law Efficiency vs. Temperature');
set(gcf, 'color', 'w');
plotfixer;

%% Part 2

T_values = [100,220,650,800] + 273; % K
lambda_range = 1:0.5:10;
P_range = 1:1:40; % atm
alpha = 0;

% Change pressure, keep lambda constant
lambda = 2;
for i=1:length(T_values)
    T = T_values(i);
    for j = 1:length(P_range)
        P = P_range(j) * P_standard;
        eta.p(i, j) = lucio(T, P, alpha, lambda);
    end
end

% Plotting for Part 2a
figure;
plot(P_range, eta.p(1, :) * 100, P_range, eta.p(2, :) * 100,...
     P_range, eta.p(3, :) * 100, P_range, eta.p(4, :) * 100);
xlabel('Pressure (Pa)');
ylabel('Efficiency (%)');
legend('80^{\circ}C', '220^{\circ}C', '650^{\circ}C', '800^{\circ}C');
title('First Law Efficiency vs. Pressure');
set(gcf, 'color', 'w');
plotfixer;

% Change lambda, keep pressure constant
P = P_standard;
for i=1:length(T_values)
    T = T_values(i);
    for j = 1:length(lambda_range)
        lambda = lambda_range(j);
        eta.lambda(i, j) = lucio(T, P, alpha, lambda);
    end
end

% Plotting for Part 2b
figure;
plot(lambda_range, eta.lambda(1, :) * 100, lambda_range, ...
     eta.lambda(2, :) * 100, lambda_range, eta.lambda(3, :) * 100,...
     lambda_range, eta.lambda(4, :) * 100)
xlabel('Lambda (-)');
ylabel('Efficiency (%)');
legend('80^{\circ}C', '220^{\circ}C', '650^{\circ}C', '800^{\circ}C');
title('First Law Efficiency vs. Lambda');
set(gcf, 'color', 'w');
plotfixer;

%% Part 3
lamda = 2; %make lambda constant again
T_values = 298:1:373;
for i = 1:length(T_values)
    T = T_values(i);
    [alpha(i), RH(i)] = relHumidity(T, lambda);
end
figure
plot(T_values,RH)
xlabel('Temperature')
ylabel('Relative Humidity')
title('Relative Humidity vs. Temperature')

    
    
    