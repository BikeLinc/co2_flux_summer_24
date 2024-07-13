clc, clear, close all

addpath(genpath('../UTILS'));

config = analysis_config();


%% Shotrcuts

deploy_7_12 = "../DATA/FLUX/07.12.2024/licor.txt";
deploy_7_12_start = datetime("7/12/2024 10:49");
deploy_7_12_end = datetime("7/12/2024 11:10");

deploy_7_10 = "../DATA/FLUX/07.10.2024/licor.txt";
deploy_7_10_start = datetime("7/10/2024 9:12");
deploy_7_10_end = datetime("7/10/2024 9:49");



%% Import Data

licor0 = IMPORTLICORFILE(deploy_7_12);
licor0 = rmmissing(licor0);
licor_idx0 = licor0.T > deploy_7_12_start & licor0.T < deploy_7_12_end;
licor0 = licor0(licor_idx0, :);

licor1 = IMPORTLICORFILE(deploy_7_10);
licor1 = rmmissing(licor1);
licor_idx1 = licor1.T > deploy_7_10_start & licor1.T < deploy_7_10_end;
licor1 = licor1(licor_idx1, :);

licor = [licor1; ];


%% smooth data
licor_sm = retime(licor,"regular", 'mean', 'TimeStep', seconds(5));
licor_sm = smoothdata(licor_sm, 'movmean', minutes(1));
licor_sm = unique(licor_sm);

%% Calculate Differences


licor_sm.DCDT = [nan;diff(licor_sm.C)];
licor_sm = rmmissing(licor_sm);

licor_sm.DCDT_mg = config.ppm_to_mg(licor_sm.DCDT);
licor_sm.C_mg = config.ppm_to_mg(licor_sm.C);

licor_sm.FLUX_mg = (licor_sm.DCDT_mg.*config.V)./config.As;


%% Expected Range

range = [0.5e-6 10e-6];% umol

range = config.mol_to_ppm(range);
range = config.ppm_to_mg(range);
disp(range);

%% Plot
figure();
plot(licor_sm.T, licor_sm.DCDT_mg,'^', 'DisplayName', 'dC/dt');
hold on;
plot(licor_sm.T, licor_sm.C_mg,'s', 'DisplayName', 'C');
legend('location','eastoutside');
xlabel("Time");
ylabel("CO_2 [mg/m^3]")
title("Measured CO_2 Emissions w/ LICOR");
grid on;

figure();
plot(licor_sm.T, licor_sm.FLUX_mg,'.', 'DisplayName', "Flux");
yline(mean(licor_sm.FLUX_mg),'r', 'DisplayName', "\mu = " + mean(licor_sm.FLUX_mg) + " mg/m^2/s")
xlabel("Time");
ylabel("Flux Rate [mg/m^2/s]");
title("Measured CO_2 Flux Rates w/ LICOR");
grid on;
legend('location','eastoutside');
yregion(range, 'FaceColor',"green", 'DisplayName', 'Expected Range');

ylim([-2 2])