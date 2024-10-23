
clc, clear, close all
%% Load Datq
load('comb_calibration_param.mat');
comb_temp_calib = daq.TA;
comb_humid_calib = daq.HA;
comb_CO2_calib = licor.C(licor.C>0);

load('amb_calibration_param.mat');
amb_temp_calib = daq.TA;
amb_humid_calib = daq.HA;
amb_CO2_calib = licor.C(licor.C>0);

load('cl_calibration_param.mat');
cl_temp_calib = daq.TA;
cl_humid_calib = daq.HA;
cl_CO2_calib = licor.C(licor.C>0);

load("measured_param.mat");
temperature_meas = daqTA;
humidity_meas = daqHA;
ambientCO2_meas = daqCA;
chamberCO2_meas = [daqCA; daqCB];

%% Fit Normal Distribution to the Data (Column 1: Combined Calibration, Column 2: Ambient Calibration, Column 3: Closed Loop Calibration, Column 4: Measured Data)
% Note: there is no ambient or chamber data for calibration, this is just one set of data
temp_mean = [mean(comb_temp_calib); mean(amb_temp_calib); mean(cl_temp_calib); mean(temperature_meas)];
temp_std = [std(comb_temp_calib); std(amb_temp_calib); std(cl_temp_calib); std(temperature_meas)];
humid_mean = [mean(comb_humid_calib); mean(amb_humid_calib); mean(cl_humid_calib); mean(humidity_meas)];
humid_std = [std(comb_humid_calib); std(amb_humid_calib); std(cl_humid_calib); std(humidity_meas)];
ambientCO2_mean = [mean(comb_CO2_calib); mean(amb_CO2_calib); mean(cl_CO2_calib); mean(ambientCO2_meas)];
ambientCO2_std = [std(comb_CO2_calib); std(amb_CO2_calib); std(cl_CO2_calib); std(ambientCO2_meas)];
chamberCO2_mean = [mean(comb_CO2_calib); mean(amb_CO2_calib); mean(cl_CO2_calib); mean(chamberCO2_meas)];
chamberCO2_std = [std(comb_CO2_calib); std(amb_CO2_calib); std(cl_CO2_calib); std(chamberCO2_meas)];


% Range of Values for Plotting (Column 1: Combined Calibration, Column 2: Ambient Calibration, Column 3: Closed Loop Calibration, Column 4: Measured Data)
x_temp = [linspace(min(comb_temp_calib) - 2, max(comb_temp_calib) + 2, 100); linspace(min(amb_temp_calib) - 2, max(amb_temp_calib) + 2, 100); linspace(min(cl_temp_calib) - 2, max(cl_temp_calib) + 2, 100); linspace(min(temperature_meas) - 2, max(temperature_meas) + 2, 100)];
x_humid = [linspace(min(comb_humid_calib) - 2, max(comb_humid_calib) + 2, 100); linspace(min(amb_humid_calib) - 2, max(amb_humid_calib) + 2, 100); linspace(min(cl_humid_calib) - 2, max(cl_humid_calib) + 2, 100); linspace(min(humidity_meas) - 2, max(humidity_meas) + 2, 100)];
x_ambientCO2 = [linspace(min(comb_CO2_calib) - 2, max(comb_CO2_calib) + 2, 100); linspace(min(amb_CO2_calib) - 2, max(amb_CO2_calib) + 2, 100); linspace(min(cl_CO2_calib) - 2, max(cl_CO2_calib) + 2, 100); linspace(min(ambientCO2_meas) - 2, max(ambientCO2_meas) + 2, 100)];
x_chamberCO2 = [linspace(min(comb_CO2_calib) - 2, max(comb_CO2_calib) + 2, 100); linspace(min(amb_CO2_calib) - 2, max(amb_CO2_calib) + 2, 100); linspace(min(cl_CO2_calib) - 2, max(cl_CO2_calib) + 2, 100); linspace(min(chamberCO2_meas) - 2, max(chamberCO2_meas) + 2, 100)];

% Calculate the Gaussian distributions (Column 1: Combined Calibration, Column 2: Ambient Calibration, Column 3: Closed Loop Calibration, Column 4: Measured Data)
pdf_temp = [normpdf(x_temp(1,:), temp_mean(1), temp_std(1)); normpdf(x_temp(2,:), temp_mean(2), temp_std(2)); normpdf(x_temp(3,:), temp_mean(3), temp_std(3)); normpdf(x_temp(4,:), temp_mean(4), temp_std(4))];
pdf_humid = [normpdf(x_humid(1,:), humid_mean(1), humid_std(1)); normpdf(x_humid(2,:), humid_mean(2), humid_std(2)); normpdf(x_humid(3,:), humid_mean(3), humid_std(3)); normpdf(x_humid(4,:), humid_mean(4), humid_std(4))];
pdf_ambientCO2 = [normpdf(x_ambientCO2(1,:), ambientCO2_mean(1), ambientCO2_std(1)); normpdf(x_ambientCO2(2,:), ambientCO2_mean(2), ambientCO2_std(2)); normpdf(x_ambientCO2(3,:), ambientCO2_mean(3), ambientCO2_std(3)); normpdf(x_ambientCO2(4,:), ambientCO2_mean(4), ambientCO2_std(4))];
pdf_chamberCO2 = [normpdf(x_chamberCO2(1,:), chamberCO2_mean(1), chamberCO2_std(1)); normpdf(x_chamberCO2(2,:), chamberCO2_mean(2), chamberCO2_std(2)); normpdf(x_chamberCO2(3,:), chamberCO2_mean(3), chamberCO2_std(3)); normpdf(x_chamberCO2(4,:), chamberCO2_mean(4), chamberCO2_std(4))];


%% Plotting

%% 4x4 Subplot with Gauysian Distributions for Calibration & Measured Data (Overlapping Combined, Ambient, Closed Loop, Measured Data) (Area Plot)
% Temp & Humidity Top Plots, Just Ambient Bottom 2 Plots
figure;

% --- Subplot 1: Temperature ---
subplot(2, 2, 1);
hold on;
grid on;
area(x_temp(1, :), pdf_temp(1, :), 'FaceColor', 'b', 'FaceAlpha', 0.2); % Combined Calibration Temperature plot in blue
area(x_temp(2, :), pdf_temp(2, :), 'FaceColor', 'r', 'FaceAlpha', 0.2); % Ambient Calibration Temperature plot in red
area(x_temp(3, :), pdf_temp(3, :), 'FaceColor', 'g', 'FaceAlpha', 0.2); % Closed Loop Calibration Temperature plot in green
area(x_temp(4, :), pdf_temp(4, :), 'FaceColor', 'm', 'FaceAlpha', 0.2); % Measured Temperature plot in magenta
xlabel('Value (%RH)');
ylabel('Probability Density');
title('Temperature Distribution');

legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Measured Data');


% --- Subplot 2: Humidity ---
subplot(2, 2, 2);
hold on;
grid on;
area(x_humid(1, :), pdf_humid(1, :), 'FaceColor', 'b', 'FaceAlpha', 0.2); % Combined Calibration Humidity plot in blue
area(x_humid(2, :), pdf_humid(2, :), 'FaceColor', 'r', 'FaceAlpha', 0.2); % Ambient Calibration Humidity plot in red
area(x_humid(3, :), pdf_humid(3, :), 'FaceColor', 'g', 'FaceAlpha', 0.2); % Closed Loop Calibration Humidity plot in green
area(x_humid(4, :), pdf_humid(4, :), 'FaceColor', 'm', 'FaceAlpha', 0.2); % Measured Humidity plot in magenta
xlabel('Value (\circC)');
ylabel('Probability Density');
title('Humidity Distribution');

legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Measured Data');

% --- Subplot 3: Ambient CO2 ---
subplot(2, 2, 3);
hold on;
grid on;
area(x_ambientCO2(1, :), pdf_ambientCO2(1, :), 'FaceColor', 'b', 'FaceAlpha', 0.2); % Combined Calibration Ambient CO2 plot in blue
area(x_ambientCO2(2, :), pdf_ambientCO2(2, :), 'FaceColor', 'r', 'FaceAlpha', 0.2); % Ambient Calibration Ambient CO2 plot in red
area(x_ambientCO2(3, :), pdf_ambientCO2(3, :), 'FaceColor', 'g', 'FaceAlpha', 0.2); % Closed Loop Calibration Ambient CO2 plot in green
area(x_ambientCO2(4, :), pdf_ambientCO2(4, :), 'FaceColor', 'm', 'FaceAlpha', 0.2); % Measured Ambient CO2 plot in magenta
xlabel('Value (ppm)');
ylabel('Probability Density');
title('Ambient CO_2 Distribution');

legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Measured Data');

% --- Subplot 4: Chamber CO2 ---
subplot(2, 2, 4);
hold on;
grid on;
area(x_chamberCO2(1, :), pdf_chamberCO2(1, :), 'FaceColor', 'b', 'FaceAlpha', 0.2); % Combined Calibration Chamber CO2 plot in blue
area(x_chamberCO2(2, :), pdf_chamberCO2(2, :), 'FaceColor', 'r', 'FaceAlpha', 0.2); % Ambient Calibration Chamber CO2 plot in red
area(x_chamberCO2(3, :), pdf_chamberCO2(3, :), 'FaceColor', 'g', 'FaceAlpha', 0.2); % Closed Loop Calibration Chamber CO2 plot in green
area(x_chamberCO2(4, :), pdf_chamberCO2(4, :), 'FaceColor', 'm', 'FaceAlpha', 0.2); % Measured Chamber CO2 plot in magenta
xlabel('Value (ppm)');
ylabel('Probability Density');
title('Chamber CO_2 Distribution');

legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Measured Data');