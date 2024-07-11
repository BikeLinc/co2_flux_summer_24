function [data, results] = CALCFLUX(data,co2_err,config, f_delivered, sp)


% convert to working units
% ppm -> mol/m^3
% lpm -> cms

co2_err = config.ppm_to_mol(co2_err);
data.C = config.ppm_to_mol(data.C);
data.CA = config.ppm_to_mol(data.CA);
data.CB = config.ppm_to_mol(data.CB);
data.Q = config.lpm_to_cms(data.Q);
config.uQ = config.lpm_to_cms(config.uQ);

% smooth dataset
data = smoothdata(data, 'movmean', config.smooth_dt);

% calculate flux
data.F = 1e6*((data.Q.*(data.CB-mean(data.CA)))./config.As);
data.F_licor = 1e6*((data.Q.*(data.C-min(data.C)))./config.As);
data.UF = 1e6*config.flux_uncert(data.CB-data.CA, config.As, co2_err, config.uQ, data.Q);

disp("Flux(licor): " + mean(data.F_licor))
disp("Flux(daq): " + mean(data.F) + " [" + mean(data.UF)+"]");
disp("Flux(delivered): " + f_delivered);


% convert to working units
% ppm -> mg/m^3
% lpm -> cms

co2_err = config.mol_to_ppm(co2_err);
data.C = config.mol_to_ppm(data.C);
data.CA = config.mol_to_ppm(data.CA);
data.CB = config.mol_to_ppm(data.CB);
co2_err = config.ppm_to_mg(co2_err);
data.C = config.ppm_to_mg(data.C);
data.CA = config.ppm_to_mg(data.CA);
data.CB = config.ppm_to_mg(data.CB);


Q = mean(data.Q);
Ca_licor = min(data.C);
Cb_licor = max(data.C);
sp = config.lpm_to_cms(sp);
f_licor = (Q*Cb_licor-(Q-sp)*Ca_licor)./config.As;
Ca_daq = min(data.CB);
Cb_daq = max(data.CB);
f_daq = (Q*Cb_daq-(Q-sp)*Ca_daq)./config.As;

f_uncert = config.flux_uncert(Cb_daq - Ca_daq, config.As, config.uQ, co2_err, Q);

f_delivered = (config.ppm_to_mg(3003)*sp)/config.As;

disp("Flux(licor): " + f_licor)
disp("Flux(daq): " + f_daq + " [" + f_uncert+"]");
disp("Flux(delivered): " + f_delivered);

results = [f_licor, f_daq, f_delivered, f_uncert];

end
