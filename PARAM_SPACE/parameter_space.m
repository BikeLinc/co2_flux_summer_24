%% parameter_space.m 
% This script is used to calculate paramater space distributions across 
% calibration datasets and field measured datasets. This file loads in
% stored .mat files containing parameterspace data that was pre-computed.

clc; clear; close all;

%% Load Files
load('comb_calibration_param.mat');
combined_temp_calib = reference.X_T;
combined_humid_calib = reference.X_H;
combined_CO2_calib = reference.Y_C;

load('amb_calibration_param.mat');
ambient_temp_calib = reference.X_T;
ambient_humid_calib = reference.X_H;
ambient_CO2_calib = reference.Y_C;

load('cl_calibration_param.mat');
closed_loop_temp_calib = reference.X_T;
closed_loop_humid_calib = reference.X_H;
closed_loop_CO2_calib = reference.Y_C;

load('measured_params.mat');
temperature_meas = [daqTA.TA(daqTA.TA>0); daqTB.TB(daqTB.TB>0)]; % remove 0 and error values.
humidity_meas = [daqHA.HA(daqHA.HA>0); daqHB.HB(daqHB.HB>0)];
ambient_CO2_meas = daqCA.CA(daqCA.CA>0);
chamber_CO2_meas = daqCB.CB(daqCB.CB>0);
ambient_CO2_meas = ambient_CO2_meas(ambient_CO2_meas~=2500); 
ambient_CO2_meas = ambient_CO2_meas(ambient_CO2_meas<5000);

%% Compute Statistics
% Temperature
temp_mean = [mean(combined_temp_calib); mean(ambient_temp_calib); mean(closed_loop_temp_calib); mean(temperature_meas)];
temp_std = [std(combined_temp_calib); std(ambient_temp_calib); std(closed_loop_temp_calib); std(temperature_meas)];

% Humidity
humid_mean = [mean(combined_humid_calib); mean(ambient_humid_calib); mean(closed_loop_humid_calib); mean(humidity_meas)];
humid_std = [std(combined_humid_calib); std(ambient_humid_calib); std(closed_loop_humid_calib); std(humidity_meas)];

% Ambient CO2
ambient_CO2_mean = [mean(combined_CO2_calib); mean(ambient_CO2_calib); mean(closed_loop_CO2_calib); mean(ambient_CO2_meas)];
ambient_CO2_std = [std(combined_CO2_calib); std(ambient_CO2_calib); std(closed_loop_CO2_calib); std(ambient_CO2_meas)];

% Chamber CO2
chamber_CO2_mean = [mean(combined_CO2_calib); mean(ambient_CO2_calib); mean(closed_loop_CO2_calib); mean(chamber_CO2_meas)];
chamber_CO2_std = [std(combined_CO2_calib); std(ambient_CO2_calib); std(closed_loop_CO2_calib); std(chamber_CO2_meas)];

%% Define Extended Plotting Ranges (Mean ± 4*Std)
% Temperature
x_temp_min = min(temp_mean - 4*temp_std);
x_temp_max = max(temp_mean + 4*temp_std);
x_temp = linspace(x_temp_min, x_temp_max, 100);

% Humidity
x_humid_min = min(humid_mean - 4*humid_std);
x_humid_max = max(humid_mean + 4*humid_std);
x_humid = linspace(x_humid_min, x_humid_max, 100);

% Ambient CO2
x_ambient_CO2_min = min(ambient_CO2_mean - 4*ambient_CO2_std);
x_ambient_CO2_max = max(ambient_CO2_mean + 4*ambient_CO2_std);
x_ambient_CO2 = linspace(x_ambient_CO2_min, x_ambient_CO2_max, 100);

% Chamber CO2
x_chamber_CO2_min = min(chamber_CO2_mean - 4*chamber_CO2_std);
x_chamber_CO2_max = max(chamber_CO2_mean + 4*chamber_CO2_std);
x_chamber_CO2 = linspace(x_chamber_CO2_min, x_chamber_CO2_max, 100);

%% Calculate Gaussian PDFs
% Temperature
pdf_temp_combined = normpdf(x_temp, temp_mean(1), temp_std(1));
pdf_temp_ambient = normpdf(x_temp, temp_mean(2), temp_std(2));
pdf_temp_closed_loop = normpdf(x_temp, temp_mean(3), temp_std(3));
pdf_temp_meas = normpdf(x_temp, temp_mean(4), temp_std(4));

% Humidity
pdf_humid_combined = normpdf(x_humid, humid_mean(1), humid_std(1));
pdf_humid_ambient = normpdf(x_humid, humid_mean(2), humid_std(2));
pdf_humid_closed_loop = normpdf(x_humid, humid_mean(3), humid_std(3));
pdf_humid_meas = normpdf(x_humid, humid_mean(4), humid_std(4));

% Ambient CO2
pdf_ambient_CO2_combined = normpdf(x_ambient_CO2, ambient_CO2_mean(1), ambient_CO2_std(1));
pdf_ambient_CO2_ambient = normpdf(x_ambient_CO2, ambient_CO2_mean(2), ambient_CO2_std(2));
pdf_ambient_CO2_closed_loop = normpdf(x_ambient_CO2, ambient_CO2_mean(3), ambient_CO2_std(3));
pdf_ambient_CO2_meas = normpdf(x_ambient_CO2, ambient_CO2_mean(4), ambient_CO2_std(4));

% Chamber CO2
pdf_chamber_CO2_combined = normpdf(x_chamber_CO2, chamber_CO2_mean(1), chamber_CO2_std(1));
pdf_chamber_CO2_ambient = normpdf(x_chamber_CO2, chamber_CO2_mean(2), chamber_CO2_std(2));
pdf_chamber_CO2_closed_loop = normpdf(x_chamber_CO2, chamber_CO2_mean(3), chamber_CO2_std(3));
pdf_chamber_CO2_meas = normpdf(x_chamber_CO2, chamber_CO2_mean(4), chamber_CO2_std(4));

%% Plot Gaussian Distributions
figure;

% Temperature
subplot(2,2,1);
hold on;
grid on;
area(x_temp, pdf_temp_combined, 'FaceColor', 'b', 'FaceAlpha', 0.2);
area(x_temp, pdf_temp_ambient, 'FaceColor', 'r', 'FaceAlpha', 0.2);
area(x_temp, pdf_temp_closed_loop, 'FaceColor', 'g', 'FaceAlpha', 0.2);
area(x_temp, pdf_temp_meas, 'FaceColor', 'm', 'FaceAlpha', 0.2);
xlabel('°C');
ylabel('Probability Density');
title('Temperature Distribution');
legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data');
xlim([0, max(x_temp)])

% Humidity
subplot(2,2,2);
hold on;
grid on;
area(x_humid, pdf_humid_combined, 'FaceColor', 'b', 'FaceAlpha', 0.2);
area(x_humid, pdf_humid_ambient, 'FaceColor', 'r', 'FaceAlpha', 0.2);
area(x_humid, pdf_humid_closed_loop, 'FaceColor', 'g', 'FaceAlpha', 0.2);
area(x_humid, pdf_humid_meas, 'FaceColor', 'm', 'FaceAlpha', 0.2);
xlabel('%RH');
ylabel('Probability Density');
title('Humidity Distribution');
legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data');
xlim([0, max(x_humid)])

% Ambient CO2
subplot(2,2,3);
hold on;
grid on;
area(x_ambient_CO2, pdf_ambient_CO2_combined, 'FaceColor', 'b', 'FaceAlpha', 0.2);
area(x_ambient_CO2, pdf_ambient_CO2_ambient, 'FaceColor', 'r', 'FaceAlpha', 0.2);
area(x_ambient_CO2, pdf_ambient_CO2_closed_loop, 'FaceColor', 'g', 'FaceAlpha', 0.2);
area(x_ambient_CO2, pdf_ambient_CO2_meas, 'FaceColor', 'm', 'FaceAlpha', 0.2);
xlabel('ppm');
ylabel('Probability Density');
title('Ambient CO_2 Distribution');
legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data');
set(gca, 'YScale', 'log');
xlim([0, max(x_ambient_CO2)])

% Chamber CO2
subplot(2,2,4);
hold on;
grid on;
area(x_chamber_CO2, pdf_chamber_CO2_combined, 'FaceColor', 'b', 'FaceAlpha', 0.2);
area(x_chamber_CO2, pdf_chamber_CO2_ambient, 'FaceColor', 'r', 'FaceAlpha', 0.2);
area(x_chamber_CO2, pdf_chamber_CO2_closed_loop, 'FaceColor', 'g', 'FaceAlpha', 0.2);
area(x_chamber_CO2, pdf_chamber_CO2_meas, 'FaceColor', 'm', 'FaceAlpha', 0.2);
xlabel('ppm');
ylabel('Probability Density');
title('Chamber CO_2 Distribution');
legend('Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data');
set(gca, 'YScale', 'log');
xlim([0, max(x_chamber_CO2)])

%% Histograms with Fitted Normal Distributions
figure;

% Define Colors and Labels
colors = {'b', 'r', 'g', 'm'};
labels = {'Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data'};

% Temperature Histogram with Fits
subplot(1,2,1);
hold on;
grid on;

% Combined Calibration
histogram(combined_temp_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{1}, 'FaceAlpha', 0.5, 'DisplayName', labels{1});
pd_comb = fitdist(combined_temp_calib, 'Normal');
x_fit = linspace(min(combined_temp_calib), max(combined_temp_calib), 100);
y_fit = pdf(pd_comb, x_fit);
plot(x_fit, y_fit, 'Color', colors{1}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

% Ambient Calibration
histogram(ambient_temp_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{2}, 'FaceAlpha', 0.5, 'DisplayName', labels{2});
pd_amb = fitdist(ambient_temp_calib, 'Normal');
y_fit_amb = pdf(pd_amb, x_fit);
plot(x_fit, y_fit_amb, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

% Closed Loop Calibration
histogram(closed_loop_temp_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{3}, 'FaceAlpha', 0.5, 'DisplayName', labels{3});
pd_cl = fitdist(closed_loop_temp_calib, 'Normal');
y_fit_cl = pdf(pd_cl, x_fit);
plot(x_fit, y_fit_cl, 'Color', colors{3}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

% Field Data
histogram(temperature_meas, 30, 'Normalization', 'pdf', 'FaceColor', colors{4}, 'FaceAlpha', 0.5, 'DisplayName', labels{4});
pd_meas = fitdist(temperature_meas, 'Normal');
y_fit_meas = pdf(pd_meas, x_fit);
plot(x_fit, y_fit_meas, 'Color', colors{4}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

xlabel('°C');
ylabel('Probability Density');
title('Temperature Histogram with Fitted Distributions');
legend('Location', 'best');

% Humidity Histogram with Fits
subplot(1,2,2);
hold on;
grid on;

% Combined Calibration
histogram(combined_humid_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{1}, 'FaceAlpha', 0.5, 'DisplayName', labels{1});
pd_comb_humid = fitdist(combined_humid_calib, 'Normal');
x_fit_humid = linspace(min(combined_humid_calib), max(combined_humid_calib), 100);
y_fit_humid = pdf(pd_comb_humid, x_fit_humid);
plot(x_fit_humid, y_fit_humid, 'Color', colors{1}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

% Ambient Calibration
histogram(ambient_humid_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{2}, 'FaceAlpha', 0.5, 'DisplayName', labels{2});
pd_amb_humid = fitdist(ambient_humid_calib, 'Normal');
y_fit_amb_humid = pdf(pd_amb_humid, x_fit_humid);
plot(x_fit_humid, y_fit_amb_humid, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

% Closed Loop Calibration
histogram(closed_loop_humid_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{3}, 'FaceAlpha', 0.5, 'DisplayName', labels{3});
pd_cl_humid = fitdist(closed_loop_humid_calib, 'Normal');
y_fit_cl_humid = pdf(pd_cl_humid, x_fit_humid);
plot(x_fit_humid, y_fit_cl_humid, 'Color', colors{3}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

% Field Data
histogram(humidity_meas, 30, 'Normalization', 'pdf', 'FaceColor', colors{4}, 'FaceAlpha', 0.5, 'DisplayName', labels{4});
pd_meas_humid = fitdist(humidity_meas, 'Normal');
y_fit_meas_humid = pdf(pd_meas_humid, x_fit_humid);
plot(x_fit_humid, y_fit_meas_humid, 'Color', colors{4}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

xlabel('%RH');
ylabel('Probability Density');
title('Humidity Histogram with Fitted Distributions');
legend('Location', 'best');

figure;

% Define Colors and Labels
colors = {'b', 'r', 'g', 'm'};
labels = {'Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data'};

% Ambient CO2 Histogram with Fits
subplot(4,1,1);
hold on;
grid on;

% Combined Calibration
histogram(combined_CO2_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{1}, 'FaceAlpha', 0.5, 'DisplayName', labels{1});
pd_comb_CO2 = fitdist(combined_CO2_calib, 'Normal');
x_fit_CO2 = linspace(min(combined_CO2_calib), max(combined_CO2_calib), 100);
y_fit_CO2 = pdf(pd_comb_CO2, x_fit_CO2);
plot(x_fit_CO2, y_fit_CO2, 'Color', colors{1}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Ambient CO_2 (ppm)');
ylabel('Probability Density');
legend('Location', 'best');
subplot(4,1,2);
hold on;
grid on;

% Ambient Calibration
histogram(ambient_CO2_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{2}, 'FaceAlpha', 0.5, 'DisplayName', labels{2});
pd_amb_CO2 = fitdist(ambient_CO2_calib, 'Normal');
y_fit_amb_CO2 = pdf(pd_amb_CO2, x_fit_CO2);
plot(x_fit_CO2, y_fit_amb_CO2, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Ambient CO_2 (ppm)');
ylabel('Log Probability Density');
legend('Location', 'best');
subplot(4,1,3);
hold on;
grid on;

% Closed Loop Calibration
histogram(closed_loop_CO2_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{3}, 'FaceAlpha', 0.5, 'DisplayName', labels{3});
pd_cl_CO2 = fitdist(closed_loop_CO2_calib, 'Normal');
y_fit_cl_CO2 = pdf(pd_cl_CO2, x_fit_CO2);
plot(x_fit_CO2, y_fit_cl_CO2, 'Color', colors{3}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Ambient CO_2 (ppm)');
ylabel('Probability Density');
legend('Location', 'best');
subplot(4,1,4);
hold on;
grid on;

% Field Data
histogram(ambient_CO2_meas, 30, 'Normalization', 'pdf', 'FaceColor', colors{4}, 'FaceAlpha', 0.5, 'DisplayName', labels{4});
pd_meas_CO2 = fitdist(ambient_CO2_meas, 'Normal');
y_fit_meas_CO2 = pdf(pd_meas_CO2, x_fit_CO2);
plot(x_fit_CO2, y_fit_meas_CO2, 'Color', colors{4}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Ambient CO_2 (ppm)');
ylabel('Probability Density');
legend('Location', 'best');

figure;

% Define Colors and Labels
colors = {'b', 'r', 'g', 'm'};
labels = {'Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data'};

% Ambient CO2 Histogram with Fits
subplot(4,1,1);
hold on;
grid on;

% Combined Calibration
histogram(combined_CO2_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{1}, 'FaceAlpha', 0.5, 'DisplayName', labels{1});
pd_comb_ChCO2 = fitdist(combined_CO2_calib, 'Normal');
x_fit_ChCO2 = linspace(min(combined_CO2_calib), max(combined_CO2_calib), 100);
y_fit_ChCO2 = pdf(pd_comb_ChCO2, x_fit_ChCO2);
plot(x_fit_ChCO2, y_fit_ChCO2, 'Color', colors{1}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Chamber CO_2 (ppm)');
ylabel('Probability Density');
legend('Location', 'best');
subplot(4,1,2);
hold on;
grid on;

% Ambient Calibration
histogram(ambient_CO2_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{2}, 'FaceAlpha', 0.5, 'DisplayName', labels{2});
pd_amb_ChCO2 = fitdist(ambient_CO2_calib, 'Normal');
y_fit_amb_ChCO2 = pdf(pd_amb_ChCO2, x_fit_ChCO2);
plot(x_fit_ChCO2, y_fit_amb_ChCO2, 'Color', colors{2}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Chamber CO_2 (ppm)');
ylabel('Log Probability Density'); 
legend('Location', 'best');
subplot(4,1,3);
hold on;
grid on;

% Closed Loop Calibration
histogram(closed_loop_CO2_calib, 30, 'Normalization', 'pdf', 'FaceColor', colors{3}, 'FaceAlpha', 0.5, 'DisplayName', labels{3});
pd_cl_ChCO2 = fitdist(closed_loop_CO2_calib, 'Normal');
y_fit_cl_ChCO2 = pdf(pd_cl_ChCO2, x_fit_ChCO2);
plot(x_fit_ChCO2, y_fit_cl_ChCO2, 'Color', colors{3}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");
xlabel('Chamber CO_2 (ppm)');
ylabel('Probability Density');
legend('Location', 'best');
subplot(4,1,4);
hold on;
grid on;

% Field Data
histogram(chamber_CO2_meas, 30, 'Normalization', 'pdf', 'FaceColor', colors{4}, 'FaceAlpha', 0.5, 'DisplayName', labels{4});
pd_meas_ChCO2 = fitdist(chamber_CO2_meas, 'Normal');
y_fit_meas_ChCO2 = pdf(pd_meas_ChCO2, x_fit_ChCO2);
plot(x_fit_ChCO2, y_fit_meas_ChCO2, 'Color', colors{4}, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', "");

xlabel('Chamber CO_2 (ppm)');
ylabel('Probability Density');
legend('Location', 'best');

% Overall Title
sgtitle('Histograms with Fitted Normal Distributions of Calibration & Field Data Parameters');

%% Chamber Violin

colors = {'b', 'r', 'g', 'm'};
labels = {'Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data'};

len1 = length(combined_CO2_calib);
len2 = length(ambient_CO2_calib);
len3 = length(closed_loop_CO2_calib);
len4 = length(chamber_CO2_meas);

% Pad vectors with nans so that can be plotted correctly.
max_len = max([len1, len2, len3, len4]);
padWithNaN = @(vec, maxLen) [vec; NaN(maxLen - length(vec), 1)];
combined_CO2_calib_padded = padWithNaN(combined_CO2_calib, max_len);
ambient_CO2_calib_padded = padWithNaN(ambient_CO2_calib, max_len);
closed_loop_CO2_calib_padded = padWithNaN(closed_loop_CO2_calib, max_len);
chamber_CO2_meas_padded = padWithNaN(chamber_CO2_meas, max_len);

data = [combined_CO2_calib_padded, ambient_CO2_calib_padded, ...
        closed_loop_CO2_calib_padded, chamber_CO2_meas_padded];

group = categorical({'Combined', 'Ambient', 'Closed Loop', 'Field Measurement'});

figure;

v = violinplot(group, data);

xlabel('Calibration Dataset');
ylabel('CO_2 Parameter Space Covered (ppm)');

grid on;

%% Ambient Violin
% Define Colors and Labels
colors = {'b', 'r', 'g', 'm'};
labels = {'Combined Calibration', 'Ambient Calibration', 'Closed Loop Calibration', 'Field Data'};


% Find the lengths of each vector
len1 = length(combined_CO2_calib);
len2 = length(ambient_CO2_calib);
len3 = length(closed_loop_CO2_calib);
len4 = length(ambient_CO2_meas);

max_len = max([len1, len2, len3, len4]);
padWithNaN = @(vec, maxLen) [vec; NaN(maxLen - length(vec), 1)];
combined_CO2_calib_padded = padWithNaN(combined_CO2_calib, max_len);
ambient_CO2_calib_padded = padWithNaN(ambient_CO2_calib, max_len);
closed_loop_CO2_calib_padded = padWithNaN(closed_loop_CO2_calib, max_len);
chamber_CO2_meas_padded = padWithNaN(ambient_CO2_meas, max_len);

data = [combined_CO2_calib_padded, ambient_CO2_calib_padded, ...
        closed_loop_CO2_calib_padded, chamber_CO2_meas_padded];

group = categorical({'Combined', 'Ambient', 'Closed Loop', 'Field Measurement'});

figure;

v = violinplot(group, data);

xlabel('Calibration Dataset');
ylabel('CO_2 Parameter Space Covered (ppm)');

grid on;

%% Combined Ambient and Chamber Violin Plot

colors = {'b', 'r', 'g', 'm', 'c', 'k'};
labels = {
    'Combined Calibration', 
    'Ambient Calibration', 
    'Closed Loop Calibration', 
    'Field Measurement Ambient', 
    'Field Measurement Chamber'
};

% Find the lengths of each vector
len1 = length(combined_CO2_calib);
len2 = length(ambient_CO2_calib);
len3 = length(closed_loop_CO2_calib);
len4 = length(ambient_CO2_meas);
len5 = length(chamber_CO2_meas);

max_len = max([len1, len2, len3, len4, len5]);
padWithNaN = @(vec, maxLen) [vec; NaN(maxLen - length(vec), 1)];
combined_CO2_calib_padded = padWithNaN(combined_CO2_calib, max_len);
ambient_CO2_calib_padded = padWithNaN(ambient_CO2_calib, max_len);
closed_loop_CO2_calib_padded = padWithNaN(closed_loop_CO2_calib, max_len);
field_CO2_meas_ambient_padded = padWithNaN(ambient_CO2_meas, max_len);
field_CO2_meas_chamber_padded = padWithNaN(chamber_CO2_meas, max_len);
data = [combined_CO2_calib_padded', 
    ambient_CO2_calib_padded', 
    closed_loop_CO2_calib_padded', 
    field_CO2_meas_ambient_padded', 
    field_CO2_meas_chamber_padded'
];

data = data';
group = categorical(labels);
figure;

v = violinplot(group, data);

xlabel('Calibration Dataset and Field Measurement');
ylabel('CO_2 Parameter Space Covered (ppm)');
title('Distribution of CO_2 Measurements Across Calibration Models and Field Data');

grid on;
set(gca, 'FontSize', 12);

