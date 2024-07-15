% This script is for running flux measurements

clc, clear, close all

addpath(genpath('../UTILS'));

cfg = config();

%% Date Range
%  Define the dates you want to look at fluxes for.
start = datetime("7/10/2024 9:12");
stop = datetime("7/10/2024 9:52");


%% Smooth Values
licorRetime = seconds(1);
licorWindow = minutes(5);
daqRetime = hours(1);
daqWindow = hours(1);

%% Collect Data
folder = sprintf('%02d',start.Month) + "." + sprintf('%02d',start.Day) + "." + start.Year;
directory = dir("../DATA/FLUX/" + folder);

licor = [];
daq = [];
for i = 3:numel(directory)
    file = directory(i);

    if file.name == "licor.txt"
        licor = IMPORTLICORFILE(file.folder + "/" + file.name);
        time = licor.T < stop - minutes(3) & licor.T > start - minutes(3);
        licor = licor(time, :);
    elseif file.name(end-3:end) == ".txt" | file.name(end-3:end) == ".TXT"
        daq = [daq; IMPORTDAQFILE(file.folder + "/" + file.name)];
        %time = daq.T < stop & daq.T > start;
        %daq = daq(time, :);
    end
end


%% Plot Raw Data

figure();
subplot(1, 2, 1);
stackedplot(licor);
subplot(1, 2, 2);
stackedplot(daq);


%% smooth data
daqRaw = daq;
licorRaw= licor;

licor = retime(licor,"regular", 'mean', 'TimeStep', licorRetime);
licor = smoothdata(licor, 'movmean', licorWindow);
licor = unique(licor);

daq = retime(daq,"regular", 'mean', 'TimeStep', daqRetime);
daq = smoothdata(daq, 'movmean', daqWindow);
daq = unique(daq);

%% Apply Calibrations

modelCA = load("../ELT CALIB/CALIB_W_LICOR/models/C-CL-linear-15-Jul-2024.mat");
modelCB = load("../ELT CALIB/CALIB_W_LICOR/models/E-CL-linear-15-Jul-2024.mat");

daq.CA_CALIB = predict(modelCA.model_1_lin, [daq.CA]);
daq.CB_CALIB = predict(modelCB.model_2_lin, [daq.CB]);

%% Plot Smoothed and Calibrated Data

figure();
subplot(1, 2, 1);
stackedplot(licor);
subplot(1, 2, 2);
stackedplot(daq);

%% Calculate Static Flux

licor.DCDT = [nan;diff(licor.C)]; % TODO:  ensure that diff returns in 
licor = rmmissing(licor);


licor.DCDT_mg = cfg.ppm_to_mg(licor.DCDT);
licor.C_mg = cfg.ppm_to_mg(licor.C);
licor.FLUX_mg = (licor.DCDT_mg.*cfg.V)./cfg.As;


%% Calculate Dynamic Flux

Q = cfg.lpm_to_cms(daq.Q/1000);
CA = cfg.ppm_to_mg(daq.CB_CALIB);
CB = cfg.ppm_to_mg(daq.CA_CALIB);

DCDT = [nan;diff(daq.CA_CALIB)];
V = cfg.V;
As = cfg.As;

daq.FLUX_mg = ((V.*DCDT) + Q.*(CB-CA))./As;
daq.FLUX_SS_mg =  (Q.*(CB-CA))./As;


%% Expected Range

range = [0.5e-6 10e-6];% umol

range = cfg.mol_to_ppm(range);
range = cfg.ppm_to_mg(range);
disp(range);

%% Plot
figure();
plot(licor.T, licor.DCDT_mg,'^', 'DisplayName', 'dC/dt');
hold on;
plot(licor.T, licor.C_mg,'s', 'DisplayName', 'C');
legend('location','eastoutside');
xlabel("Time");
ylabel("CO_2 [mg/m^3]")
title("Measured CO_2 Emissions w/ LICOR");
grid on;

figure();
hold on;
plot(licor.T, licor.FLUX_mg,'.', 'DisplayName', "LICOR");
plot(daq.T, daq.FLUX_SS_mg,'.', 'DisplayName', "DAQ");
yline(mean(licor.FLUX_mg),'r', 'DisplayName', "LICOR \mu = " + mean(licor.FLUX_mg) + " mg/m^2/s")
yline(mean(rmmissing(daq.FLUX_SS_mg)),'', 'DisplayName', "DAQ \mu = " + mean(rmmissing(daq.FLUX_SS_mg)) + " mg/m^2/s")
xlabel("Time");
ylabel("Flux Rate [mg/m^2/s]");
title("Measured CO_2 Flux Rates w/ LICOR");
grid on;
legend('location','eastoutside');
yregion(range, 'FaceColor',"green", 'DisplayName', 'Expected Range');

%ylim([-2 2])

%%

figure()
hold on;
plot(daq.T, daq.CA, 'DisplayName', 'CA');
plot(daq.T, daq.CA_CALIB,'--', 'DisplayName', 'CA CALIB');
plot(daq.T, daq.CB, 'DisplayName', 'CB');
plot(daq.T, daq.CB_CALIB,'--', 'DisplayName', 'CB CALIB');
legend();