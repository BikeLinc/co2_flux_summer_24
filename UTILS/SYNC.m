function [data] = SYNC(daq, licor, map, config, dataset)
%SYNC This function synchronizes datasets into a single timetable and
%applies and offset to daq datasets that is manually determined.

    % get offset from map table
    sp_idx = table2array(map(:,1)) == double(dataset);
    map = map(sp_idx,:);
    daqoffset = min(table2array(map(sp_idx,8)));

    % daq dataset was set to epoch, resetting shifting by start of flux tests
    if dataset == "5.21"
        daq.T = timeofday(daq.T) + daqoffset - minutes(69);
    elseif dataset == "5.22"
        daq.T = timeofday(daq.T) + daqoffset - minutes(34);
    elseif dataset == "5.29"
        daq.T = daq.T + hours(1) + minutes(42);
    end
    
    % downsample datasets
    daq = retime(daq, 'regular', 'mean', 'TimeStep', config.sample_dt);
    licor = retime(licor, 'regular', 'mean', 'TimeStep', config.sample_dt);
    
    % synchronize datasets, and cleanup dataset  
    data = synchronize(daq, licor);
    best_corr = 0;
    opt_lag = -inf;
    for lag = -height(data):height(data)
        if lag > 0
            shifted_CB = [nan(lag, 1); data.CB(1:end-lag)];
        elseif lag < 0
            shifted_CB = [data.CB(-lag+1:end); nan(-lag, 1)];
        else
            shifted_CB = data.CB;
        end
        
        % Calculate correlation, ignoring NaNs
        valid_idx = ~isnan(data.C) & ~isnan(shifted_CB);
        if sum(valid_idx) > 100
            current_corr = corr(data.C(valid_idx), shifted_CB(valid_idx));
            % Update best correlation and lag
            if current_corr > best_corr
                best_corr = current_corr;
                opt_lag = lag;
            end
        end
        
        
    end
     fields = ["CA", "CB", "TA", "TB", "HA", "HB","Q"];
    for field = fields
        data.(field) = SHIFTDATA(data.(field), opt_lag);
    end
    data = rmmissing(data);
end

