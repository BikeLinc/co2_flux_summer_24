%% lab_analysis.m
%
%   This file is for interpereting delivered fluxes and validating the CO2
%   flux chamber in summer 2024. We use the licor, measuring the chamber
%   exhaust concentration and the flux chambers measurement, calculate the
%   flux, then report the plots.
%
%   Lincoln Scheer
%   Jun 23 2024
%
%   Updated: 10/13/2024
%   - Added linear regression to correct for bias in low-cost system
%   - Added plots to show corrected fluxes

clc, clear, close all

addpath(genpath('../UTILS'));

config = analysis_config();

%% Load Data & Process
% Import data from daq (low-cost) and licor (LI-7810) and synchronize datasets with set-points delivered.
%   - Import data
%   - Synchronize data
%   - Apply calibrations
%   - Calculate fluxes
%   - Store results
%   - Report results
%   - Store fluxes

fluxes = [];
timeseries_fluxes = [];
flux_count = 0;

for dataset = config.datasets
    
    % import dataset
    [daq_data, licor_data] = IMPORTDATA(config.path, dataset);
    set_point_map = readtable(config.map_path);
    set_point_map = set_point_map(table2array(set_point_map(:,1)) == double(dataset),:);

    % synchronize dataset
    data = SYNC(daq_data, licor_data, set_point_map, config, dataset);

    % apply calibrations
    [data, co2_error] = CALIBRATE(data, dataset);

    for set_point_index = 1:1:height(set_point_map)
        
        % get set point timestamps
        licor_start_time = set_point_map{set_point_index, 10};
        licor_end_time = set_point_map{set_point_index, 11};

        % get delivered fluxes
        flux_delivered = set_point_map{set_point_index, 13};

        % select corresponding setpoint
        flow_rate = set_point_map{set_point_index, 12};
        data_set_point_index = data.T < licor_end_time  & data.T > licor_start_time;
        data_set_point = data(data_set_point_index, :);        

        % calculate fluxes
        [data_set_point, results] = CALCFLUX(data_set_point, co2_error, config, flux_delivered, flow_rate);

        % store dataset
        writetimetable(data_set_point, "data/analysis/" + dataset + "no" + set_point_index +"results.csv");

        % report results
        %disp(results);
        disp("-")

        % store fluxes
        fluxes = [fluxes; results];

        % store timeseries fluxes
        flux_count = flux_count + 1;
        timeseries_fluxes(flux_count).data = data_set_point;
        timeseries_fluxes(flux_count).results = results;
        
    end

end



%% Plot Results

% sort flux responses by delivered flux
fluxes = sortrows(fluxes);

% extract fluxes
delivered = fluxes(:,3);
low_cost = fluxes(:,2);
low_cost_err = fluxes(:,4);
li_7810 = fluxes(:,1);

% convert units to mg/m^2/min
delivered = delivered * 60;
low_cost = low_cost * 60;
low_cost_err = low_cost_err * 60;
li_7810 = li_7810 * 60;

% generate linear regression to correct for bias. Map low-cost to LI-7810
% and correct for bias
lm = fitlm(low_cost, li_7810);
bias = lm.Coefficients.Estimate(1);
slope = lm.Coefficients.Estimate(2);

% plot 1:1 fit and corrected fluxes
figure();
hold on;
plot(delivered,delivered, '--', 'DisplayName', '1:1 Fit');
errorbar(delivered, low_cost*slope + bias, low_cost_err, 'x-', 'DisplayName', "Low-Cost System w/ Correction");
errorbar(delivered, low_cost, low_cost_err, '.-', 'DisplayName', "Low-Cost System");
plot(delivered, li_7810, 'o-.','DisplayName', "LI-7810");
xlabel("Delivered CO_2 Flux [mg/m^2/min]")
ylabel("Measured CO_2 Flux [mg/m^2/min]")
grid on;
legend('Location', 'southeast');

% print equation
disp("Correction Equation: y = " + slope + "x + " + bias);
disp("R^2: " + lm.Rsquared.Ordinary);