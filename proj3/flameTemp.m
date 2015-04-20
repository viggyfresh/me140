function [ T_a ] = flameTemp(phi, type, hf_mol, MM)

T_a = 25 + 273.15;
dT = 0.1;

if (strcmp(type, 'JetA'))
    lhs = phi * hf_mol.JetA;
    rhs = 12.3 * phi * hf_mol.CO2 + 11.1 * phi * hf_mol.H2O;
    
    while rhs < lhs
        T_a = T_a + dT;
        rhs = rhs + 12.3 * phi * (sp_heats(T_a, 'CO2') / 1000 * MM.CO2 * dT);
        rhs = rhs + 11.1 * phi * (sp_heats(T_a, 'H2O') / 1000 * MM.H2O * dT);
        coeff = (2 * 12.3 + 11.1) / 2;
        rhs = rhs + (coeff * (79 / 21) * sp_heats(T_a, 'N2') / 1000 * MM.N2 * dT);
        rhs = rhs + (coeff * (1 - phi) * sp_heats(T_a, 'O2') / 1000 * MM.O2 * dT);
    end
elseif (strcmp(type, 'Dodecane'))
    lhs = phi * hf_mol.Dodecane;
    rhs = 12 * phi * hf_mol.CO2 + 13 * phi * hf_mol.H2O;
    
    while rhs < lhs
        T_a = T_a + dT;
        rhs = rhs + 12 * phi * (sp_heats(T_a, 'CO2') / 1000 * MM.CO2 * dT);
        rhs = rhs + 13 * phi * (sp_heats(T_a, 'H2O') / 1000 * MM.H2O * dT);
        coeff = (2 * 12 + 13) / 2;
        rhs = rhs + (coeff * (79 / 21) * sp_heats(T_a, 'N2') / 1000 * MM.N2 * dT);
        rhs = rhs + (coeff * (1 - phi) * sp_heats(T_a, 'O2') / 1000 * MM.O2 * dT);
    end
end