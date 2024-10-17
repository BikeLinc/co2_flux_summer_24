%% Standard functions for manuscript calculations
%
% Lincoln Scheer
% 8/21/2024
%

run var.m


f.ppm_to_mol = @(ppm) (ppm*P)/(1e6*8.3145*T);       % ppm to mol/m^3
f.mol_to_ppm = @(mol) (1e6*8.3145*mol*T)/P;         % mol/m^3 to ppm
f.ppm_to_mg = @(ppm) (ppm*P*44010)/(1e6*8.3145*T);  % ppm to mg/m^3
f.mg_to_ppm = @(mol) (1e6*8.3145*mol*T)/(P*44010);  % mg/m^3 to ppm
f.lpm_to_cms = @(lpm) lpm/60000;                    % liters per min to m^3 per min
f.cms_to_lpm = @(cms) cms*60000;                    % m^3 per min to liters per min