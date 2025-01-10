clc; clear; close all;

addpath(genpath('../UTILS'));
config = analysis_config();
calibration_methods = {'linear', 'ANN'};
calibration_results = struct();
fluxes_all_methods = struct();

for method_idx = 1:length(calibration_methods)
    method = calibration_methods{method_idx};
    fprintf('Starting analysis with %s calibration...\n', method);
    
    fluxes = [];
    timeseries_fluxes = [];
    flux_count = 0;
    
    for dataset = config.datasets
        [daq_data, licor_data] = IMPORTDATA(config.path, dataset);
        set_point_map = readtable(config.map_path);
        set_point_map = set_point_map(double(set_point_map{:,1}) == double(dataset), :);
        
        data = SYNC(daq_data, licor_data, set_point_map, config, dataset);
        [data_calibrated, uncalib_data, co2_error] = CALIBRATE2(data, dataset, method);

        fig = figure();
        hold on;
        grid on;
        plot(data_calibrated.C, 'b-', 'DisplayName', "Reference") 
        plot(data_calibrated.CB,'g-','DisplayName', "CO_2 Chamber [Calibration Applied]")
        plot(uncalib_data.CB,'g.', 'DisplayName', "CO_2 Chamber [Uncalibrated]")
        legend()
        xlabel("Time (min)")
        ylabel("CO_2 (ppm)")

        fig_filename = sprintf('lab_flux_calibration_timeseries_%s_%s.png', method, dataset);
        saveas(fig, fig_filename)
        close all;

        for set_point_index = 1:height(set_point_map)
            licor_start_time = set_point_map{set_point_index, 10};
            licor_end_time = set_point_map{set_point_index, 11};
            flux_delivered = set_point_map{set_point_index, 13};
            flow_rate = set_point_map{set_point_index, 12};
            data_set_point_idx = data_calibrated.T < licor_end_time & data_calibrated.T > licor_start_time;
            data_set_point = data_calibrated(data_set_point_idx, :);        
    
            [data_set_point, results] = CALCFLUX(data_set_point, co2_error, config, flux_delivered, flow_rate);
    
            output_filename = sprintf('data/analysis/%s_no%d_results_%s.csv', dataset, set_point_index, method);
            writetimetable(data_set_point, output_filename);
    
            disp(results);
            disp('-');
    
            fluxes = [fluxes; results];
            flux_count = flux_count + 1;
            timeseries_fluxes(flux_count).data = data_set_point;
            timeseries_fluxes(flux_count).results = results;
        end
    end
    
    fluxes = sortrows(fluxes);
    
    % Extract fluxes
    delivered = fluxes(:,3);
    low_cost = fluxes(:,2);
    low_cost_err = fluxes(:,4);
    li_7810 = fluxes(:,1);
    
    % Convert units to mg/m^2/min
    delivered = delivered * 60;
    low_cost = low_cost * 60;
    low_cost_err = low_cost_err * 60;
    li_7810 = li_7810 * 60;
    
    % Generate linear regression to correct for bias: Map low-cost to LI-7810
    lm = fitlm(low_cost, li_7810);
    bias = lm.Coefficients.Estimate(1);
    slope = lm.Coefficients.Estimate(2);
    
    calibration_results.(method).slope = slope;
    calibration_results.(method).bias = bias;
    calibration_results.(method).Rsquared = lm.Rsquared.Ordinary;
    
    fluxes_all_methods.(method).delivered = delivered;
    fluxes_all_methods.(method).low_cost = low_cost;
    fluxes_all_methods.(method).low_cost_err = low_cost_err;
    fluxes_all_methods.(method).li_7810 = li_7810;
    
    figure('Name', sprintf('Calibration Results - %s', method));
    hold on;
    grid on;
    
    plot(delivered, delivered, 'k-', 'DisplayName', '1:1 Fit');
    errorbar(delivered, low_cost * slope + bias, low_cost_err, '.--', 'DisplayName', 'Low-Cost System w/ Correction');
    plot(delivered, li_7810, '.-.', 'DisplayName', 'LI-7810');
    
    xlabel('Delivered CO_2 Flux [mg/m^2/min]');
    ylabel('Measured CO_2 Flux [mg/m^2/min]');
    title(sprintf('Flux Calibration using %s Method', method));
    legend('Location', 'southeast');
    
    % Display Calibration Equation and R²
    disp(['Calibration Method: ', method]);
    disp(['Correction Equation: y = ' num2str(slope) 'x + ' num2str(bias)]);
    disp(['R²: ', num2str(lm.Rsquared.Ordinary)]);
    
    % Save Figure
    fig_filename = sprintf('lab_flux_calibration_%s.png', method);
    saveas(gcf, fig_filename);
    close(gcf);
end

figure('Name', 'Combined Calibration Results - Linear vs ANN');
hold on;
grid on;

plot(delivered, delivered, 'k-', 'DisplayName', '1:1 Fit', 'LineWidth', 2);


x_fit = linspace(min(delivered), max(delivered), 100);

% Linear Regression
linear_slope = calibration_results.linear.slope;
linear_bias = calibration_results.linear.bias;
linear_R2 = calibration_results.linear.Rsquared;

y_linear_corrected = fluxes_all_methods.linear.low_cost * linear_slope + linear_bias;

errorbar(fluxes_all_methods.linear.delivered, y_linear_corrected, fluxes_all_methods.linear.low_cost_err, ...
         'rx', 'DisplayName', sprintf('Fluxes (R²=%.2f)', linear_R2), 'MarkerSize', 5, 'LineWidth', 1.5);

lm_lin = fitlm(delivered, y_linear_corrected);
plot(delivered, delivered*lm_lin.Coefficients.Estimate(2)+lm_lin.Coefficients.Estimate(1),'r-.', ...
     'DisplayName', "Fluxes Best Fit", 'LineWidth',2)

xlabel('Delivered CO_2 Flux [mg/m^2/min]');
ylabel('Measured CO_2 Flux [mg/m^2/min]');
legend('Location', 'northwest');

pbaspect([1 1 1])
saveas(gcf, 'lab_flux_calibration_LR.fig');
close(gcf);

figure;
hold on;
grid on;

plot(delivered, delivered, 'k-', 'DisplayName', '1:1 Fit', 'LineWidth', 2);

ANN_slope = calibration_results.ANN.slope;
ANN_bias = calibration_results.ANN.bias;
ANN_R2 = calibration_results.ANN.Rsquared;

y_ANN_corrected = fluxes_all_methods.ANN.low_cost * ANN_slope + ANN_bias;

errorbar(fluxes_all_methods.ANN.delivered, y_ANN_corrected, fluxes_all_methods.ANN.low_cost_err, ...
         'bo', 'DisplayName', sprintf('Fluxes (R²=%.2f)', ANN_R2), 'MarkerSize', 5, 'LineWidth', 1.5);

lm_ann = fitlm(delivered, y_ANN_corrected);
plot(delivered, delivered*lm_ann.Coefficients.Estimate(2)+lm_ann.Coefficients.Estimate(1),'b-.', ...
     'DisplayName', "Fluxes Best Fit", 'LineWidth',2)

xlabel('Delivered CO_2 Flux [mg/m^2/min]');
ylabel('Measured CO_2 Flux [mg/m^2/min]');
legend('Location', 'northwest');

pbaspect([1 1 1])
saveas(gcf, 'lab_flux_calibration_ANN.fig');
close(gcf);

fprintf('Analysis completed for all calibration methods.\n');