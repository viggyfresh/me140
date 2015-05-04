function [ cp, cv, gamma, R ] = sp_heats( temp , type )

% Molecular masses - grams
MM.O2 = 32;
MM.N2 = 28.02;
MM.C = 12.01;
MM.H = 1.008;
MM.H2 = 2 * MM.H;
MM.H2O = 18.016;
MM.CO2 = MM.C + MM.O2;
MM.air = 28.97;

if strcmp(type,'CO2')
    mm = MM.CO2;
    a = 22.26;
    b = 5.981*10^-2;
    c = -3.501*10^-5;
    d = 7.469*10^-9;
    
elseif strcmp(type,'H2O_vap')
    mm = MM.H2O;
    a = 32.24;
    b = 0.1923*10^-2;
    c = 1.055*10^-5;
    d = -3.595*10^-9;
    
elseif strcmp(type, 'H2O_liq')
    mm = 1;
    a = 4.18;
    b = 0;
    c = 0;
    d = 0;
    
elseif strcmp(type,'N2')
    mm = MM.N2;
    a = 28.9;
    b = -0.1571*10^-2;
    c = 0.8081 * 10^-5;
    d = -2.873*10^-9;
    
elseif strcmp(type,'O2')
    mm = MM.O2;
    a = 25.48;
    b = 1.520*10^-2;
    c = -0.7155*10^-5;
    d = 1.312*10^-9;
    
elseif strcmp(type,'H2')
    mm = MM.H2;
    a = 29.11;
    b = -0.1916e-2;
    c = 0.4003e-5;
    d = -.8704e-9;

elseif strcmp(type,'air')
    mm = MM.air;
    a = 28.11;
    b = 0.1967*10^-2;
    c = 0.4802*10^-5;
    d = -1.966*10^-9;
end
    
% J / kg everywhere!
R = 8.314462 / mm * 1000;
p = [d c b a];
cp = (polyval(p,temp)) / mm * 1000;
cv = cp - R;
gamma = cp ./ cv;

end