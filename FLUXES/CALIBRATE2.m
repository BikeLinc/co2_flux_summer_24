function [data, uncalib_data, co2_err] = CALIBRATE2(data, dataset, method)
    % CALIBRATE2 applies calibration to the data using the specified method
    % and retains the uncalibrated data with '_UNCALIB' appended to the variable names.
    %
    % Inputs:
    %   data     - The data table containing raw measurements.
    %   dataset  - Identifier for the dataset (e.g., "5.29").
    %   method   - Calibration method: 'linear' or 'ANN'.
    %
    % Outputs:
    %   data     - The calibrated data table with uncalibrated data appended.
    %   co2_err  - The maximum CO2 error from calibration.

    % Store uncalibrated data with '_UNCALIB' appended to variable names
    uncalib_data = data;

    % Load calibration models
    load('calib.mat', '*');

    switch lower(method)
        case 'linear'
            % Apply linear regression calibrations
            if dataset == "5.29"
                data.CB = predict(lin_regb, data.CB);
                data.CA = predict(lin_rega, data.CA);
            else
                data.CB = predict(lin_rega, data.CB);
                data.CA = predict(lin_regb, data.CA);
            end

            % Apply manual linear regressions
            data.C = data.C .* 0.9883 + 24.6914;
            data.Q = data.Q * 1.227 + 0.0143;

            % Set CO2 error based on linear regression RMSE
            co2_err = rmse(data.CB,data.C);

        case 'ann'
            % Apply ANN regression calibrations
            if dataset == "5.29"
                data.CB = ann_regb(data.CB')';
                data.CA = ann_rega(data.CA')';
            else
                data.CB = ann_rega(data.CB')';
                data.CA = ann_regb(data.CA')';
            end

            % Apply manual linear regressions if needed
            data.C = data.C .* 0.9883 + 24.6914;
            data.Q = data.Q * 1.227 + 0.0143;

            % Set CO2 error based on ANN regression RMSE

            
            
            co2_err = rmse(data.CB,data.C);
        otherwise
            error('Unknown calibration method. Choose "linear" or "ANN".');
    end
end
