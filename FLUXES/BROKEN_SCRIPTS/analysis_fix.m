%% Flux Test Data Analysis - 5/27/2024
% This script is for analyzing flux datasets with the CO2 Flux Chamber
% designed by the Wildfire CO2 Flux Team at Fort Lewis College.
%
% This script is desinged to work with a LICOR LI-7810 reference sensor.
%

clc, clear, close all

addpath(genpath('../UTILS'));

%% Script Settings

% Chamber Surface Area
As = 0.064844702; % [m^2]
V  = 0.015220183; % [m^3]

% Environment Assumptions
P = 78082.9; % [Pa]
T = 293.15; % [K]
MW = 44010; % [mg/mol]

% Flow Meas. Uncertainty
uQ = 0.33; % [lpm]

% MFC Uncertainty (%RD)
uQ_mfcperc = 2;

% Preprocessingw
sample_dt = seconds(5);     % retiming applied to entire dataset
smooth_dt = minutes(2);   % retiming applied to per set-point dataset

% Choose Dataset (File Location)
dataset = "5.29";

%% Helper Functions

% Unit Conversion Functions
ppm_to_mgm3 = @(ppm) (ppm*P*MW)/(1e6*8.3145*T);  % ppm to mol/m^3
mgm3_to_ppm = @(mol) (1e6*8.3145*mol*T)/(P*MW);    % mol/m^3 to ppm
lpm_to_cms = @(lpm) lpm/60000;              % liters per min to m^3 per min
cms_to_lpm = @(cms) cms*60000;              % m^3 per min to liters per min

% Uncertainty in Flux Measurement Function
flux_uncert = @(dC, As, uQ, uC, Q) sqrt( (dC.*uQ./As).^2 + ((2.*Q.*uC)./As).^2 );


%% Import Datasets

% import licor dataset
licor = IMPORTLICORFILE('data/'+dataset+'/licor.data');

% import daq dataset
daq = IMPORTDAQFILE('data/'+dataset+'/daq.txt');
daq.Q = daq.Q/1000; % flow-meter from ccm to lpm

% import delivered flux dataset
map = readtable("data/mapping.csv");

% depending on dataset folder selected, map delivered fluxes
sp_idx = table2array(map(:,1)) == double(dataset);
map = map(sp_idx,:);
daqoffset = min(table2array(map(:,8)));

% daq dataset was set to epoch, resetting shifting by start of flux tests
if dataset == '5.21'
daq.T = timeofday(daq.T) + daqoffset - minutes(69);
elseif dataset == '5.22'
daq.T = timeofday(daq.T) + daqoffset - minutes(34);
elseif dataset == '5.29'
daq.T = daq.T + hours(1) + minutes(42);
end

% sychronize both datasets
sync_data = synchronize(licor, daq, 'regular', 'mean','TimeStep', sample_dt);

clear licor daq
%% Correct Timestamps
%if dataset ~= '5.29'
%sync_data = correctTimestamps(sync_data);
%end
%% Calibrate Sensors

[sync_data, ca_rmse, cb_rmse] = applyCalibrations(sync_data, 0, dataset);


%% Seperate Datasets into Set Points
metrics_sp = zeros(height(map), 4);
data_sp = cell(height(map),1);

for sp_idx = 1:height(map)
    % get sp timestamps
    tStart = map{sp_idx, 8};
    tStartLicor = map{sp_idx, 10};
    tEnd = map{sp_idx, 9};
    tEndLicor = map{sp_idx, 11};
    
    % select setpoint
    data_idx = sync_data.T < tEndLicor  & sync_data.T > tStartLicor;
    data = sync_data(data_idx, :);
    
    % apply retime average
    %data = retime(data, 'regular', 'mean', 'TimeStep', smooth_dt);
    %data = rmmissing(data);
    data = smoothdata(data, 'movmean', smooth_dt);


    % calculate theoretical steady state
    [tss_data_l, tss_l, Cchmb_l] = CO2CHAMBERTSS(1, ppm_to_mgm3(mean(data.CA_CALIB)), map{sp_idx, 13}*1e-6, As, lpm_to_cms(mean(data.Q_CALIB)), 1, V, ppm_to_mgm3(5));
    tss_data_l.Time = seconds(tss_data_l.TIME) + data.T(1);
    tss_data_l.CO2= mgm3_to_ppm(tss_data_l.CO2);

    % calculate floor dataset, because we are looking for offsets
    data.C_FLOOR = data.C_CALIB - min(data.C_CALIB);
    data.CB_FLOOR = data.CB_CALIB - min(data.CB_CALIB);
    data.CA_FLOOR = data.CA_CALIB - min(data.CA_CALIB);

    % calculate steady state indices for ELT Chamber and LICOR
    cb_ss = findSteadyStateIndices(data, "CB_CALIB", cb_rmse);
    c_ss = findSteadyStateIndices(data, "C_CALIB", cb_rmse);

    % calculate flux



    data.F = ((lpm_to_cms(data.Q_CALIB).*ppm_to_mgm3(data.CB_CALIB-data.CA_CALIB))./As).*1e6;

    data.F_LICOR = ((lpm_to_cms(data.Q_CALIB).*ppm_to_mgm3(data.C -data.C(1)))./As).*1e6;
        
    data.UF = flux_uncert(ppm_to_mgm3(data.CB-data.CA), As, lpm_to_cms(uQ), ppm_to_mgm3(cb_rmse+ca_rmse/2), lpm_to_cms(data.Q))*1e6;
    data.UF_LICOR = flux_uncert(ppm_to_mgm3(data.CB-data.CA), As, lpm_to_cms(uQ), ppm_to_mgm3(1.5), lpm_to_cms(data.Q))*1e6;

    f_mean_ss = mean(data.F(end));
    uf_mean_ss = mean(data.UF(end));
    f_licor_mean_ss = mean(data.F_LICOR(end));
    uf_licor_mean_ss = mean(data.UF_LICOR(end));

    metrics_sp(sp_idx, :) = [f_mean_ss, uf_mean_ss, f_licor_mean_ss, uf_licor_mean_ss];
    
    figure();
   
    % plot raw set-point dataset
    subplot(2, 3, 1)
    hold on
    plot(data.T, data.C, 'b.-', 'DisplayName', 'LICOR Reference')
    plot(data.T, data.CB, 'g.-', 'DisplayName', 'DAQ Chamber')
    plot(data.T, data.CA, 'r.-', 'DisplayName', 'DAQ Ambient')
    ylabel('CO_2 [ppm]');
    legend();
    title(["Raw Data - Delivering " + map{sp_idx, 13} + " μmol/m^2/s", "[DATASET " + dataset + "]"]);
    xlabel('Time');
    hold off
    
    % plot calibrated set-point dataset
    subplot(2, 3, 2)
    hold on
    plot(data.T, data.C_CALIB, 'b.-', 'DisplayName', 'Corrected LICOR Reference')
    plot(data.T, data.CB_CALIB, 'g.-', 'DisplayName', 'Corrected DAQ Chamber')
    plot(data.T, data.CA_CALIB, 'r.-', 'DisplayName', 'Corrected DAQ Ambient')
    plot(tss_data_l.Time, tss_data_l.CO2, 'k-', 'DisplayName', 'Expected Chamber Concentration', 'LineWidth', 2)
    ylabel('CO_2 [ppm]');
    legend();
    title(["Calibrated Data - Delivering " + map{sp_idx, 13} + " μmol/m^2/s", "[DATASET " + dataset + "]"]);
    xlabel('Time');
    hold off

    % plot floored dataset
    subplot(2, 3, 3)
    hold on;
    plot(data.T, data.CB_FLOOR, 'g.-', 'DisplayName', 'DAQ CB');
    plot(data.T, data.CA_FLOOR, 'r.-', 'DisplayName', 'DAQ CA');
    plot(data.T, data.C_FLOOR, 'c.-', 'DisplayName', 'LICOR')
    plot(tss_data_l.Time, tss_data_l.CO2-min(tss_data_l.CO2), 'k-', 'DisplayName', 'Expected Chamber Concentration', 'LineWidth', 2)
    ylabel('CO_2 [ppm]');
    legend();
    title(["Floored Data - Delivering " + map{sp_idx, 13} + " μmol/m^2/s", "[DATASET " + dataset + "]"]);
    xlabel('Time');
    hold off;
    
    % plot steady state indices
    subplot(2, 3, 4)
    hold on;
    plot(data.T(c_ss), data.CB_CALIB(c_ss), 'gd', 'DisplayName', 'Steady-State DAQ Chamber');
    plot(data.T(cb_ss), data.CA_CALIB(cb_ss), 'rd', 'DisplayName', 'Steady-State DAQ Ambient');
    plot(data.T(cb_ss), data.C_CALIB(cb_ss), 'cd', 'DisplayName', 'Steady-State LICOR')
    plot(data.T, data.CB_CALIB, 'g.-', 'DisplayName', 'DAQ Chamber');
    plot(data.T, data.CA_CALIB, 'r.-', 'DisplayName', 'DAQ Chamber');
    plot(data.T, data.C_CALIB, 'c.-', 'DisplayName', 'LICOR')
    plot(tss_data_l.Time, tss_data_l.CO2, 'k-', 'DisplayName', 'Expected Chamber Concentration', 'LineWidth', 2)
    ylabel('CO_2 [ppm]');
    legend();
    title(["Steady State Indices - Delivering " + map{sp_idx, 13} + " μmol/m^2/s", "[DATASET " + dataset + "]"]);
    xlabel('Time');
    hold off
    % plot fluxes
    subplot(2, 3, 5)
    hold on;
    errorbar(data.T(cb_ss), data.F(cb_ss), data.UF(cb_ss), 'g.-', 'DisplayName', "SS Flux, " + f_mean_ss + " μmol/m^2/s");
    errorbar(data.T(c_ss), data.F_LICOR(c_ss), data.UF_LICOR(c_ss), 'r.-', 'DisplayName', "SS Flux LICOR, " + f_licor_mean_ss + " μmol/m^2/s");
    yline(map{sp_idx, 13}, 'k--', 'DisplayName', "Delivered Flux")
    ylabel('CO_2 Flux μmol/m^2/s');
    legend();
    title(["Steady State Flux Results - Delivering " + map{sp_idx, 13} + " μmol/m^2/s", "[DATASET " + dataset + "]"]);
    xlabel('Time');
    hold off;

    % plot fluxes
    subplot(2, 3, 6)
    hold on;
    plot(data.T, data.F, 'g.-', 'DisplayName', "Flux");
    plot(data.T, data.F_LICOR, 'r.-', 'DisplayName', "Flux LICOR");
    yline(map{sp_idx, 13}, 'k--', 'DisplayName', "Delivered Flux")
    ylabel('CO_2 Flux μmol/m^2/s');
    legend();
    title(["Flux Results - Delivering " + map{sp_idx, 13} + " μmol/m^2/s", "[DATASET " + dataset + "]"]);
    xlabel('Time');
    hold off;
    disp(data.C_CALIB(1));
    disp(data.C_CALIB(end));
    
end
%%
figure()

subplot(1, 2, 1)
hold on
errorbar(map{:, 13}, metrics_sp(:,1), metrics_sp(:, 2), metrics_sp(:, 2), map{:, 14}, map{:, 14},'.-', 'DisplayName', 'ELT Chamber Flux')
plot(map{:, 13}, map{:, 13}, 'r--', 'DisplayName', '1:1 Fit')
ylabel('Measured CO_2 Flux [μmol/m^2/s]');
xlabel('Delivered CO_2 Flux [μmol/m^2/s]');
legend()
title("DAQ Flux Results")
grid on
xlim("padded")
ylim("padded")
subplot(1, 2, 2)

hold on
errorbar(map{:, 13}, metrics_sp(:,3), metrics_sp(:, 4), metrics_sp(:, 4), map{:, 14}, map{:, 14},'.-', 'DisplayName', 'LICOR Flux')
plot(map{:, 13}, map{:, 13}, 'r--', 'DisplayName', '1:1 Fit')
ylabel('Measured CO_2 Flux [μmol/m^2/s]');
xlabel('Delivered CO_2 Flux [μmol/m^2/s]');
legend()
title("LICOR Flux Results")
grid on
xlim("padded")
ylim("padded")
sgtitle("Comparing Measured to Delivered Fluxes")

figure()

plot(map{:, 13}, metrics_sp(:,1)./metrics_sp(:, 2), '.-','DisplayName', 'Signal to Noise Ration')
xlabel('Delivered CO_2 Flux [μmol/m^2/s]');
ylabel('Signal to Noise Measured CO_2 Flux')
title("Signal to Noise Ratio of Measured Fluxes vs. Delivered Fluxes")
grid on
xlim("padded")
ylim("padded")

%% Functions for Data Analysis

function corrected_data = correctTimestamps(sync_data)
    %% Correct Timestamps (Only If DAQ Lost RTC Data)

    % if daq has timestamp issue, do offset with start of test sequence


    % Find Cross-Corelation Lag & Offset Datasets

    % synchronize dataset, resample at n-seconds per datapoint
    
    sync_data_c = sync_data.C;
    sync_data_cb = sync_data.CB;
    
    % computer lag-based corelation coefficient
    opt_lag = 0;
    opt_corr = -inf;
    
    max_lag = min(850, floor(height(sync_data) / 8));
    
    for lag = -max_lag:max_lag
        
        % offset dataset
        if lag > 0
            sync_data_cb_shifted = [nan(lag, 1); sync_data_cb(1:end-lag)];
        elseif lag < 0
            sync_data_cb_shifted = [sync_data_cb(-lag+1:end); nan(-lag, 1)];
        else
            sync_data_cb_shifted = sync_data_cb;
        end
        
        % calculation corellation
        non_nan_indices = ~isnan(sync_data_c) & ~isnan(sync_data_cb_shifted);
        current_corr = corr(sync_data_c(non_nan_indices), sync_data_cb_shifted(non_nan_indices));
    
        % update best corelation
        if current_corr > opt_corr
            opt_corr = current_corr;
            opt_lag = lag;
        end
    end

    fields_to_shift = {'CB', 'CA', 'TA', 'TB', 'HA', 'HB', 'Q'};
    corrected_data = sync_data;
    for i = 1:length(fields_to_shift)
        field = fields_to_shift{i};
        corrected_data.(field) = shiftData(sync_data.(field), opt_lag);
    end
    corrected_data = rmmissing(corrected_data);
    
end

function shift_data = shiftData(data, lag)

    if lag > 0
        shift_data = [nan(lag, 1); data(1:end-lag)];
    elseif lag < 0
        shift_data = [data(-lag+1:end); nan(-lag, 1)];
    else
        shift_data = data;
    end

end

function [calibrated_data, ca_rmse, cb_rmse] = applyCalibrations(sync_data, verbose, dataset)

    load('calib.mat', '*');

    if dataset == '5.29'
        sync_data.CB_CALIB = predict(lin_regb, sync_data.CB);
        sync_data.CA_CALIB = predict(lin_rega, sync_data.CA);
    else
        sync_data.CB_CALIB = predict(lin_rega, sync_data.CB);
        sync_data.CA_CALIB = predict(lin_regb, sync_data.CA);
    end
   
    sync_data.C_CALIB = sync_data.C.*0.9883+24.6914;
    sync_data.Q_CALIB = (sync_data.Q*1.227+0.0143);
    % apply ANN regressions, instead of linear
    %corr_data.CB_CALIB = ann_regb(corr_data.CB')';
    %corr_data.CA_CALIB = ann_rega(corr_data.CA')';
    
    if verbose
        figure();
        hold on;
        plot(sync_data.T, sync_data.CB_CALIB, 'g.-', 'DisplayName', 'Corrected DAQ Chamber');
        plot(sync_data.T, sync_data.CA_CALIB, 'r.-', 'DisplayName', 'Corrected DAQ Ambient');
        plot(sync_data.T, sync_data.C_CALIB, 'b.-', 'DisplayName', 'Corrected LICOR')
        plot(sync_data.T, sync_data.CB, 'c.-', 'DisplayName', 'DAQ Chamber')
        plot(sync_data.T, sync_data.CA, 'm.-', 'DisplayName', 'DAQ Ambient')
        plot(sync_data.T, sync_data.C, 'k.-', 'DisplayName', 'LICOR')
        ylabel('CO_2 [ppm]');
        legend();
        title(["Calibrations Applied to Flux Dataset" "[DATASET]"]);
        xlabel('Time');
        hold off;
    end

    calibrated_data = sync_data;
    ca_rmse = lin_rega.RMSE;
    cb_rmse = lin_regb.RMSE;
end

function tss_idx = findSteadyStateIndices(data, var, thresh)

    tss = data.T(1);
    for i = 2:length(data.(var))
        dC = abs(data.(var)(i) - data.(var)(i-1));
        if(dC > thresh)
            tss = data.T(i);
        end
    end

    tss_idx = data.T > tss;

    if tss_idx(end) == 0
        tss_idx(end) = 1;
    end

end