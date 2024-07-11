%%  Dilution Analysis
%
%   This script is for analyzing dilution test data collected on 4 April
%   2024 in the M&I lab.
%
%   CO2 Flux Project
%
%   5 April 2024
%   Lincoln Scheer
%
%

clc, clear, close all


%% Importing Data

% Importing DAQ Dataset
daq = IMPORTDAQ2("data/April 16 2024/daq.csv");
daq = rmmissing(daq);
daq = table2timetable(daq);

% Importing LICOR Dataset
licor = IMPORTLICOR("data/April 16 2024/licor.txt");
licor.T = datetime( licor.DATE + timeofday(licor.TIME) , 'Format', 'default');
licor.DATE = [];
licor.TIME = [];
licor = rmmissing(licor);
licor = table2timetable(licor);

% Importing Dilution Mapping Dataset
test = IMPORTTEST("data/April 16 2024/mapping.csv");
test_start = min(test.START);
test_end = max(test.END);

% Confine LICOR Dataset to Test Range
idx_licor = licor.T < test_end & licor.T > test_start;
licor = licor(idx_licor, :);

% Confine LICOR Dataset to Test Range
idx_daq = daq.T < test_end & daq.T > test_start;
daq = daq(idx_daq, :);




%% Preprocessing

% Apply a 10-Seconds Retime w/ Linear Interpolation
retime_dt = seconds(30);
licor = retime(licor, 'regular', 'linear', 'TimeStep', retime_dt);
daq = retime(daq, 'regular', 'linear', 'TimeStep', retime_dt);

clear idx_licor retime_dt

%% Mapping Dilutions( LICOR and DAQ ), Adding TARGET Column to Tables

licor_map = zeros(height(licor),1);

% LICOR Mapping
licor.TARGET = zeros(height(licor),1);
for i = 1:height(test)
    range = licor.T > test.START(i) & licor.T < test.END(i);
    licor.TARGET(range) = ones(sum(range),1)*test.TARGET(i);
    licor_map(range) = ones(sum(range),1)*test.TARGET(i);
end

licor.TARGET(licor.TARGET == 1) = NaN;
licor.TARGET(licor.TARGET == 0) = NaN;
licor_map_tt = timetable(licor.T, licor_map, 'VariableNames', "TARGET");
licor = rmmissing(licor);
licor_map_tt = rmmissing(licor_map_tt);

% DAQ Mapping
daq_map = zeros(height(daq),1);
daq.TARGET = zeros(height(daq),1);
for i = 1:height(test) 
    range = daq.T > test.START(i) & daq.T < test.END(i);
    daq.TARGET(range) = ones(sum(range),1)*test.TARGET(i);
    daq_map(range) = ones(sum(range),1)*test.TARGET(i);
end
daq.TARGET(daq.TARGET == 1) = NaN;
daq.TARGET(daq.TARGET == 0) = NaN;

daq_map_tt = timetable(daq.T, daq_map, 'VariableNames', "TARGET");

daq = rmmissing(daq);
daq_map_tt = rmmissing(daq_map_tt);


%% Regression Finding

licor_lm = fitlm(licor.C, licor.TARGET, 'PredictorVars', 'CO2 Licor', 'ResponseVar', 'CO2 Target');
daq_lma = fitlm(daq.CA, daq.TARGET, 'PredictorVars', 'CO2 DAQ A', 'ResponseVar', 'CO2 Target');
daq_lmb = fitlm(daq.CB, daq.TARGET, 'PredictorVars', 'CO2 DAQ B', 'ResponseVar', 'CO2 Target');

models = {licor_lm, daq_lma, daq_lmb};

%% Plotting Data

fig1 = figure();
hold on
plot(licor.T, licor.C);
plot(licor.T, licor.TARGET);
plot(daq.T, daq.CA);
plot(daq.T, daq.CB);
legend(["LICOR", "TARGET", "DAQ A", "DAQ B"])
xlabel("Timestamp")
ylabel("CO_2 [ppm]")
title("Timeseries of Dilution Test Variables")
grid on

%% Plotting Regressions

for i = 1:length(models)
    figure()
    subplot(1, 2, 1)
    plot(models{:, i})
    title(models{:, i}.VariableNames(1) + " Fit")
    subplot(1, 2, 2)
    plotResiduals(models{i}, 'fitted')
    title(models{:, i}.VariableNames(1) + " Fitted Residuals")
end

%%


for i = 1:length(models)
    
    disp(models{:, i}.VariableNames(1) + " RMSE: " + models{:, i}.RMSE)
    
end

daqa_std = std(daq.CA)
daqb_std = std(daq.CB)
licor_std = std(licor.C)

%%

plot(daq_lma)

%% Plotting Error Bars

% DAQ A
sem = mean(daq.CA)/sqrt(height(daq));
eem = 2*sem;
error_bars = ones(size(daq.T))*eem;
figure

shadedErrorBar(datenum(daq.T) , daq.CA, error_bars, 'lineprops', '-r');
%dateaxis('x', 14, daq.T(1))
xlabel("Timestamp")
ylabel("CO_2 [ppm]")
title("ELT A CO2 Measurement and 95% CI (Shaded)")
grid on
hold off

figure

errorbar(datenum(daq.T) , daq.CA, error_bars);
%dateaxis('x', 14, daq.T(1))
xlabel("Timestamp")
ylabel("CO_2 [ppm]")
title("ELT A CO2 Measurement and 95% CI (Shaded)")
grid on
hold off

%% Error

datesa = [min(test.START), min(test.END)];
idxa = daq.T < datesa(2) & daq.T > datesa(1);
daqa = daq(idxa, :);

datesb = [max(test.START), max(test.END)];
idxb = daq.T < datesb(2) & daq.T > datesb(1);
daqb = daq(idxb, :);

% DAQ A
sema = mean(daqa.CA)/sqrt(height(daqa));
eema = 2*sema;
error_barsa = ones(size(daqa.T))*eema;

semb = mean(daqb.CA)/sqrt(height(daqb));
eemb = 2*semb;
error_barsb = ones(size(daqb.T))*eemb;

meana = mean(daqa.CA);
meanb = mean(daqb.CA);

figure()

errorbar([height(daqa), height(daqb)], [meana, meanb], [eema, eemb], "o", 'LineWidth', 2)
xlabel("Target Dilevered CO_2[ppm]")
ylabel("Mean Measured Sensor CO_2[ppm]")
title("Comparison of Mean Measured Sensor CO_2 At Set-Point")