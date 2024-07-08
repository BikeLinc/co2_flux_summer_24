%% This Dataset is For Closed Loop Sensor Calibration
%
% Lincoln Scheer
% 6/16/2024
%

clc, clear, close all

%% Import Data

% import elt sensor dataset
daq = IMPORTDAQFILE("data/daq_data_7_5.csv");
daq = rmmissing(daq);

% import licor reference instrument dataset
licor = IMPORTLICORFILE("data/licor_data_7_5.txt");
licor = rmmissing(licor);

% split per-sensor data
sensor_data = {[daq(:,[2,3,4])], [daq(:,[5,6,7])]};


%%
% remove rows that show the ELTs throwing errors
elt_errors = [-999,500, 2815, 64537, 231753, 65535, 2500, 2559];

% remove values that are outside of expected dynamic range.


%% plot raw sensor data
figure();
hold on;
plot(daq.T, daq.CA, '.','DisplayName', 'ELT CO_2 A');
plot(daq.T, daq.CB, '.','DisplayName', 'ELT CO_2 B');
plot(daq.T, daq.TA, '.','DisplayName', 'TEMP A');
plot(daq.T, daq.HA, '.','DisplayName', 'HUMID A');
plot(daq.T, daq.TA, '.','DisplayName', 'TEMP B');
plot(daq.T, daq.HA, '.','DisplayName', 'HUMID B');
plot(licor.T, licor.C, '.','DisplayName', 'LICOR CO_2');
legend();
title("Closed Loop Calibration - Raw Sensor Data");

%% smooth data over 5 mins

smooth_dt = minutes(1);
retime_dt = seconds(5);

% retime datasets
daqa = retime(daqa,"regular", 'mean', 'TimeStep', retime_dt);
daqb = retime(daqb,"regular", 'mean', 'TimeStep', retime_dt);
licor = retime(licor,"regular", 'mean', 'TimeStep', retime_dt);

% apply moving mean
%daqa = rmoutliers(daqa, 'percentile', [20 80]);
%daqb = rmoutliers(daqb, 'percentile', [20 80]);
%licor = rmoutliers(licor, 'percentile', [20 80]);

% apply moving mean
daqa = smoothdata(daqa, 'movmean', smooth_dt);
daqb = smoothdata(daqb, 'movmean', smooth_dt);
licor = smoothdata(licor, 'movmean', smooth_dt);



dataa = synchronize(daqa, licor);
datab = synchronize(daqb, licor);

dataa = rmmissing(dataa);
datab = rmmissing(datab);

% plot smoothed data

figure();
hold on;
plot(dataa.T, dataa.CA, 'DisplayName', 'ELT CO_2 A');
plot(datab.T, datab.CB, 'DisplayName', 'ELT CO_2 B');
plot(licor.T, licor.C, 'DisplayName', 'LICOR CO_2');
legend();
title("Closed Loop Calibration - Hourly Smoothed Sensor Data");

%% generate regressions

CA_Model = fitlm([dataa.CA, dataa.TA, dataa.HA], dataa.C);
CB_Model = fitlm([datab.CB, datab.TB, datab.HB], datab.C);


% plot regressions
text_size = 30;

fig = figure();
hold on

plot(CA_Model.Variables.y, CA_Model.Variables.y, '--', 'LineWidth', 5, 'MarkerSize', 20);
plot(CA_Model.Fitted, CA_Model.Variables.y, '.', 'LineWidth', 5, 'MarkerSize', 20);


xlabel("ELT A CO_2 [ppm]",'Interpreter','tex');
ylabel("LICOR CO_2 [ppm]", 'Interpreter','tex');
title('Chamber Sensor Calibration','Interpreter','tex');
legend(["1:1 Fit","Fitted CO_2 Dataset"], 'Interpreter', 'tex');
fontsize(fig, 30, 'points')
fontname('Times New Roman')

txt = "RMSE: " + round(CA_Model.RMSE,3) + " ppm\newlineR^2: " + round(CA_Model.Rsquared.Ordinary,3) + "\newliney="+round(table2array(CA_Model.Coefficients(1,2)),1)+ "x_1+"+round(table2array(CA_Model.Coefficients(1,3)),1)+ "x_2+"+round(table2array(CA_Model.Coefficients(1,4)),1)+ "x_3+"+round(table2array(CA_Model.Coefficients(1,1)),1);

disp("Model A");
disp(txt);

%text(min(xlim)+5, max(ylim)-30,  txt,'Interpreter','tex', 'FontSize', 40, 'FontName', 'Times New Roman');


%



fig2 = figure();
hold on

plot(CB_Model.Variables.y, CB_Model.Variables.y, '--', 'LineWidth', 5, 'MarkerSize', 20);
plot(CB_Model.Fitted, CB_Model.Variables.y, '.', 'LineWidth', 5, 'MarkerSize', 20);


xlabel("ELT A CO_2 [ppm]",'Interpreter','tex');
ylabel("LICOR CO_2 [ppm]", 'Interpreter','tex');
title('Chamber Sensor Calibration','Interpreter','tex');
legend(["1:1 Fit","Fitted CO_2 Dataset"], 'Interpreter', 'tex');
fontsize(fig2, 30, 'points')
fontname('Times New Roman')

txt = "RMSE: " + round(CB_Model.RMSE,3) + " ppm\newlineR^2: " + round(CB_Model.Rsquared.Ordinary,3) + "\newliney="+round(table2array(CB_Model.Coefficients(1,2)),1)+ "x_1+"+round(table2array(CB_Model.Coefficients(1,3)),1)+ "x_2+"+round(table2array(CB_Model.Coefficients(1,4)),1)+ "x_3+"+round(table2array(CB_Model.Coefficients(1,1)),1);
disp("Model B");
disp(txt);

%text(min(xlim)+5, max(ylim)-30,  txt,'Interpreter','tex', 'FontSize', 40, 'FontName', 'Times New Roman');