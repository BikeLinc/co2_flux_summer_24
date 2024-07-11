%% modular_analysis.m
%
%
%   This file is for interpereting delivered fluxes and validating the CO2
%   flux chamber in summer 2024.
%
%   Lincoln Scheer
%   Jun 23 2024
%
%

clc, clear, close all

%% Calculate Fluxes

% load configuration variables
config = analysis_config();

% calculate per-dataset fluxes, each dataset consists of a collection os
% set-points that we delivered.
fluxes = [];
for dataset = config.datasets
    
    % import dataset
    [daq, licor] = IMPORTDATA(config.path, dataset);
    map = readtable(config.map_path);
    map = map(table2array(map(:,1)) == double(dataset),:);

    % synchronize dataset
    data = SYNC(daq, licor, map, config, dataset);

    % apply calibrations
    [data, co2_err] = CALIBRATE(data, dataset);

    

    for sp_idx = 1:1:height(map)
        
        % get sp timestamps
        tStartLicor = map{sp_idx, 10};
        tEndLicor = map{sp_idx, 11};

        % get delivered fluxes
        f_delivered = map{sp_idx, 13};

        % select setpoint
        sp = map{sp_idx, 12};
        data_sp_idx = data.T < tEndLicor  & data.T > tStartLicor;
        data_sp = data(data_sp_idx, :);        

        figure()
        hold on
        plot(data_sp.T, data_sp.C,'.');
        plot(data_sp.T, data_sp.CB,'.');
        plot(data.T, data.C,'o');
        plot(data.T, data.CB,'o');


        % calculate fluxes
        [data_sp, results] = CALCFLUX(data_sp, co2_err, config, f_delivered, sp);

        % store dataset
        writetimetable(data_sp, "data/analysis/" + dataset + "no" + sp_idx +"results.csv");

        % report results
        %disp(results);
        disp("-")

        % store fluxes
        fluxes = [fluxes; results];

    end

end

%% Plot

fluxes = sortrows(fluxes);

fig = figure();
hold on;
plot(fluxes(:,3),fluxes(:,3), '--', 'DisplayName', '1:1 Fit', 'LineWidth', 5);
errorbar(fluxes(:,3), fluxes(:,2), fluxes(:,4), 'x-', 'DisplayName', "Low-Cost System", 'LineWidth', 5, 'MarkerSize', 10, 'CapSize', 5);
plot(fluxes(:,3), fluxes(:,1), 'o-.','DisplayName', "LI-7810", 'LineWidth', 5, 'MarkerSize', 20);
title("Laboratory Flux Measurements");
xlabel("Delivered CO_2 Flux [mg/m^2/s]")
ylabel("Measured CO_2 Flux [mg/m^2/s]")
legend()
fontsize(fig, 50, 'points')
fontname('Times New Roman')