function [data, co2_err] = CALIBRATE(data, dataset)
    load('calib.mat', '*');

    if dataset == "5.29"
        data.CB = predict(lin_regb, data.CB);
        data.CA = predict(lin_rega, data.CA);
    else
        data.CB = predict(lin_rega, data.CB);
        data.CA = predict(lin_regb, data.CA);
    end
   
    % manual linear regressions
    data.C = data.C.*0.9883+24.6914;
    data.Q = data.Q*1.227+0.0143;

    % apply ANN regressions, instead of linear
    %corr_data.CB_CALIB = ann_regb(corr_data.CB')';
    %corr_data.CA_CALIB = ann_rega(corr_data.CA')';

    co2_err = max(lin_rega.RMSE,lin_regb.RMSE);
end