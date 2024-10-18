%% closed_loop.m
% This script is for generating calibration models for the closed-loop test
% between the ELT S300 sensor and LICOR LI-7810 reference instrument. This
% script does not support multiple calibration files currently.
%
% Lincoln Scheer
% 7/8/2024
%
% This script has not been optimized and can take a while to run, sit back.
%

clc, clear, close all

addpath(genpath('../../UTILS'));

%% Load Data


%% Sensor Data (ELT E = 1, ELT C = 2)
% Closed Loop Data
cl_files = ["../../DATA/CALIB/07.08.2024/daq_7_8_24.csv", "../../DATA/CALIB/07.08.2024/daq_7_9_24.csv"];
cl_times = ["7/8/2024  9:45:00" "7/9/2024  7:00:00"];
cl_daq = [];
for i = 1:length(cl_files)
    daq = IMPORTDAQFILE(cl_files(i));
    idx = daq.T > datetime(cl_times(1)) & daq.T < datetime(cl_times(2));
    cl_daq = [cl_daq; daq(idx, :)];
end

% Ambient Data
amb_files = ["../../DATA/CALIB/07.09.2024/24_07_09.TXT", "../../DATA/CALIB/07.09.2024/24_07_10.TXT"];
amb_times = ["7/9/2024  7:30:00" "7/10/2024  6:40:00"];
amb_daq = [];
for i = 1:length(amb_files)
    daq = IMPORTDAQFILE(amb_files(i));
    idx = daq.T > datetime(amb_times(1)) & daq.T < datetime(amb_times(2));
    amb_daq = [amb_daq; daq(idx, :)];
end

% Combine Datasets
daq = [cl_daq; amb_daq];
%daq = [amb_daq];
daq = rmmissing(daq);

% from elt sensor dataset, grab per-sensor dataset
% sensor 1 & 2 (CA & CB)

%% Reference Data
% Closed Loop Data
cl_licor = IMPORTLICORFILE("../../DATA/CALIB/07.08.2024/licor.txt");
cl_licor_idx = cl_licor.T > datetime(cl_times(1)) & cl_licor.T < datetime(cl_times(2));
cl_licor = cl_licor(cl_licor_idx, :);

% Ambient Data
amb_licor = IMPORTLICORFILE("../../DATA/CALIB/07.09.2024/licor_7_10_24.txt");
amb_licor_idx = amb_licor.T > datetime(amb_times(1)) & amb_licor.T < datetime(amb_times(2));
amb_licor = amb_licor(amb_licor_idx, :);

% Combine Datasets
licor = [cl_licor; amb_licor];
%licor = [amb_licor];
licor.T = licor.T + minutes(3); % LICOR timestamp is 3 minutes behind
licor = rmmissing(licor);

%% Combine Datasets into Per-Sensor Timetables
sensors = {[daq(:,[2,3,4])], [daq(:,[5,6,7])]};

%% Remove Errors
errorVals = [-999, 500, 2815, 64537, 231753, 65535, 2500, 2559];

% Remove known error values from each dataset
for index = 1:2
    sensor = sensors{1,index};
    sensor = timetable2table(sensor);
    errorMask = (ismember(table2array(sensor(:,2)), errorVals));
    sensor = sensor(~errorMask, :);
    sensor = table2timetable(sensor);
    sensors{1,index} = sensor;
end

%% Plot Raw Data
figure();
hold on; grid on;
plot(daq.T, daq.CA,'DisplayName', 'ELT CO_2 A', 'Marker', 'diamond', 'LineStyle', 'none');
plot(daq.T, daq.CB,'DisplayName', 'ELT CO_2 B', 'Marker', '^', 'LineStyle', 'none');
plot(licor.T, licor.C,'DisplayName', 'LICOR CO_2', 'Marker', 'square', 'LineStyle', 'none');
xlabel("Time")
ylabel("CO_2 [ppm]")
title("Closed Loop Calibration - Raw Sensor Data");
legend('location','eastoutside')

%% Retime, Smooth, and Remove Outliers
% section settings
smooth_dt = minutes(15);
retime_dt = seconds(10);
outlier_bounds = [2, 98];
outlier_remove = true;

% smooth and retime sensor datasets
for index = 1:2
    sensor = sensors{1,index};
    sensor = retime(sensor,"regular", 'mean', 'TimeStep', retime_dt);
    sensor = smoothdata(sensor, 'movmean', smooth_dt);
    if (outlier_remove)
        sensor = rmoutliers(sensor, 'percentile', outlier_bounds);
    end
    sensors{1,index} = sensor;
end

% smooth and retime reference dataset
licor = retime(licor,"regular", 'mean', 'TimeStep', retime_dt);
licor = smoothdata(licor, 'movmean', smooth_dt);
if (outlier_remove)
    licor = rmoutliers(licor, 'percentile', outlier_bounds);
end

%% Synchronize Datasets

% sync sensors with reference
for index = 1:2
    sensor = sensors{1,index};
    sensor = synchronize(sensor, licor);
    sensor = rmmissing(sensor);

    sensor = renamevars(sensor, [1, 2, 3, 4], ["X_C", "X_T", "X_H", "Y_C"]);

    sensors{1,index} = sensor;
end

%% Interaction Terms

for index = 1:2
    sensor = sensors{1,index};


    sensor.X_CDT = sensor.X_C ./ sensor.X_T;
    sensor.X_CDH = sensor.X_C ./ sensor.X_H;
    sensor.X_LNC = log(sensor.X_C);
    sensor.X_CDLNC = sensor.X_C ./ log(sensor.X_C);

    sensors{1,index} = sensor;
end

%% On-The-Fly Timestamp Fix

% look for timestamp lags, and automatically fix based on cross corelation
% of lagged datasets
for index = 1:2
    sensor = sensors{1,index};
    best_corr = 0;
    opt_lag = -inf;

    for lag = -height(sensor):height(sensor)
        if lag > 0
            shifted_C = [nan(lag, 1); sensor.(1)(1:end-lag)];
        elseif lag < 0
            shifted_C = [sensor.(1)(-lag+1:end); nan(-lag, 1)];
        else
            shifted_C = sensor.(1);
        end
        
        % calculate correlation, ignoring NaNs
        valid_idx = ~isnan(sensor.Y_C) & ~isnan(shifted_C);
        if sum(valid_idx) > 100
            current_corr = corr(sensor.Y_C(valid_idx), shifted_C(valid_idx));
            % update best correlation and lag
            if current_corr > best_corr
                best_corr = current_corr;
                opt_lag = lag;
            end
        end
    end

    % apply best lag and shift dataset
    fields = 1:3; 
    for field = fields
        tmp_shft = SHIFTDATA(sensor.(field), opt_lag);
        sensor.(field) = tmp_shft;
    end
    sensor = rmmissing(sensor);
end

%% Plot Smooth Data
figure();
sgtitle("Closed Loop Calibration - Processed Sensor Data");
for index = 1:2
    subplot(1,2,index);
    sensor = sensors{1,index};
    hold on; grid on;
    plot(sensor, 1, 'DisplayName', 'ELT', 'Marker', 'diamond', 'LineStyle', 'none');
    plot(sensor, 4, 'DisplayName','LICOR', 'Marker', 'square', 'LineStyle', 'none');
    xlabel("Time");
    ylabel("CO_2 [ppm]")
    legend();
    sensors{1,index} = sensor;
end

%% Generate Calibrations

predictors = ["X_C", "X_T", "X_CDT"]%, "X_CDLNC"];
target = "Y_C";

models = cell(2,2);
for index = 1:2
    sensor = sensors{1,index};
    sensor = timetable2table(sensor);

    % partition data
    cv = cvpartition(size(sensor,1), 'HoldOut', 0.3);
    trainX = table2array(sensor(~cv.test, predictors)); %2:4 for all predictors
    trainY = table2array(sensor(~cv.test, target));
    testX = table2array(sensor(cv.test, predictors)); %2:4 for all predictors
    testY = table2array(sensor(cv.test, target));

    % linear regression
    models{1,index} = fitlm(trainX, trainY);

    % network regression (with tolerance of 0.1ppm)
    models{2,index} = feedforwardnet([16, 16]);
    models{2,index}.trainParam.min_grad = 1;
    models{2,index}.trainParam.epochs = 1000;
    models{2,index}.trainParam.max_fail = 10;
    models{2,index}.trainParam.goal = 10;
    models{2,index} = train(models{2,index}, trainX', trainY');

    % calculate metrics
    lin_pred = predict(models{1,index}, testX);
    net_pred = models{2,index}(testX')';
    lin_r2 = 1 - ((sum((lin_pred - testY).^2))/(sum(((testY - mean(testY)).^2))));
    net_r2 = 1 - ((sum((net_pred - testY).^2))/(sum(((testY - mean(testY)).^2))));
    lin_rmse = models{1,index}.RMSE;
    net_rmse = sqrt(mean((net_pred - testY).^2));

    sensor = table2timetable(sensor);
    sensors{1,index} = sensor;
end



%%
fig_resid = figure()
fprintf("Closed Loop Calibration - Calibration Results\n")

cnt = 0;
for index = 1:2
    sensor = sensors{1,index};
    sensor = timetable2table(sensor);

    % plot metrics - (w/ Dr. Casey suggestions)

    cnt = cnt + 1;
    subplot(2,2,cnt)
    hold on; grid on;
    plot(testY, testY, 'r-.', 'DisplayName', "1:1 Fit");
    plot(testY, lin_pred ,'gs', 'DisplayName', "Linear Model");
    plot(testY, net_pred,'mo', 'DisplayName', "Neural Model");
    legend('location','best');
    xlabel("  CO_2 [ppm]")
    ylabel("Model CO_2 [ppm]")
    title(".")
    
    cnt = cnt + 1;
    subplot(2,2,cnt)
    hold on; grid on;
    yline(0, 'r-.', 'DisplayName', "1:1 Fit")
    plot(testY, lin_pred-testY,'gs', 'DisplayName', "Linear Model Response Residuals");
    plot(testY, net_pred-testY,'mo', 'DisplayName', "Neural Model Response Residuals");
    legend('location','best');
    xlabel("True CO_2 [ppm]")
    ylabel("True CO_2 Residuals [ppm]")
    title(".")
 
    % print metrics
    fprintf("Sensor %d\n", index);
    fprintf("\tLinear\n")
    fprintf("\t\tRMSE:\t%0.2f\n", lin_rmse);
    fprintf("\t\tR^2:\t%0.5f\n", lin_r2);
    fprintf("\tNetwork\n")
    fprintf("\t\tRMSE:\t%0.2f\n", net_rmse);
    fprintf("\t\tR^2:\t%0.5f\n", net_r2);

    sensor = table2timetable(sensor);
    sensors{1,index} = sensor;
    cnt = 2;
end


%% Export Models

sensor_id = "C";
training_id = "CL";
model_1_lin = models{1,1};
model_1_net = models{2,1};
save("models/"+sensor_id+"-"+training_id+"-linear-"+datestr(datetime("today")),"model_1_lin");
save("models/"+sensor_id+"-"+training_id+"-net-"+datestr(datetime("today")),"model_1_net");

sensor_id = "E";
training_id = "CL";
model_2_lin = models{1,2};
model_2_net = models{2,2};
save("models/"+sensor_id+"-"+training_id+"-linear-"+datestr(datetime("today")),"model_2_lin");
save("models/"+sensor_id+"-"+training_id+"-net-"+datestr(datetime("today")),"model_2_net");

%% Add non training data performance

% use ambient data to test model performance
% prepare ambient data inputs
amb_sensors = {[amb_daq(:,[2,3,4])], [amb_daq(:,[5,6,7])]};
amb_licor.T = amb_licor.T + minutes(3);

fig_perf = figure();


for index = 1:2
    subplot(2,2,index);
    sensor = sensors{1,index};
    hold on; grid on;
    plot(sensor, 1, 'DisplayName', 'Raw NDIR Sensor', 'Marker', '.', 'LineStyle', 'none');
    plot(sensor, 4, 'DisplayName','LI-7810', 'Marker', '.', 'LineStyle', 'none');
    title("A")
    xlabel("Time");
    ylabel("CO_2 [ppm]")
    legend('Location', 'best');
    sensors{1,index} = sensor;
end

for index = 1:2
    sensor = amb_sensors{1,index};
    sensor = retime(sensor,"regular", 'mean', 'TimeStep', retime_dt);
    sensor = smoothdata(sensor, 'movmean', smooth_dt);
    sensor = rmoutliers(sensor, 'percentile', outlier_bounds);
    sensor = timetable2table(sensor);
    sensor = renamevars(sensor, [1, 2, 3, 4], ["T", "X_C", "X_T", "X_H"]);

    sensor.X_CDT = sensor.X_C ./ sensor.X_T;
    sensor.X_CDH = sensor.X_C ./ sensor.X_H;
    sensor.X_LNC = log(sensor.X_C);
    sensor.X_CDLNC = sensor.X_C ./ log(sensor.X_C);

    prediction = predict(models{1,index}, table2array(sensor(:, predictors)));
    sensor.Y_C_LM = prediction;
    prediction = models{2,index}(table2array(sensor(:, predictors))')';
    sensor.Y_C_NN = prediction;

    % plot results
    subplot(2, 2, index + 2);
    hold on; grid on;
    plot(sensor.T, sensor.X_C, '.', 'DisplayName', 'Raw NDIR Sensor');
    plot(sensor.T, sensor.Y_C_LM, '.','DisplayName', 'w/ Linear Model Applied ');
    plot(sensor.T, sensor.Y_C_NN, '.', 'DisplayName', 'w/ FFN Model Applied');
    yline(426, 'r-.', 'DisplayName', 'NASA July 2024 Ambient CO_2');
    xlabel("Time")
    ylabel("CO_2 [ppm]")
    title("Ambient Calibration - Sensor " + index)
    legend('Location', 'best');

    amb_sensors{1,index} = sensor;
end

%% Generate Calibration 1:1 Plots for Validation Data
