function [ hf_mol, hf_kg ] = heatOfFormation()


LHV_kg.JetA = 42800 * 10^3; % J / kg
LHV_kg.Dodecane = 44560 * 10^3; % J / kg

LHV_mol.JetA = LHV_kg.JetA * 170.145 / 1000; % J / mol
LHV_mol.Dodecane = LHV_kg.Dodecane * 170.38 / 1000; % J / mol

hf_mol.H2O = -241820; % vapor, J / mol 
hf_mol.CO2 = -393520; % J / mol
hf_mol.JetA = 12.3 * hf_mol.CO2 + 11.1 * hf_mol.H2O + LHV_mol.JetA;
hf_mol.Dodecane = 13 * hf_mol.H2O + 12 * hf_mol.CO2 + LHV_mol.Dodecane;

hf_kg.H2O = hf_mol.H2O / 18.02 * 1000;
hf_kg.CO2 = hf_mol.CO2 / 44.01 * 1000;
hf_kg.JetA = hf_mol.JetA / 170.145 * 1000;
hf_kg.Dodecane = hf_mol.Dodecane / 170.38 * 1000;