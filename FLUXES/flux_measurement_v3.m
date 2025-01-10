clc; clear; close all;

fprintf("Running:\tflux_measurement.m\n");
addpath(genpath('../UTILS'));
cfg = config();

dateStart = datetime(2024, 7, 8);
dateEnd = datetime(2024, 7, 31);
startTime = timeofday(datetime("00:59", 'InputFormat', 'HH:mm')) + dateStart;
endTime = timeofday(datetime("11:59", 'InputFormat', 'HH:mm')) + dateEnd;

processStatic = false;
processDynamic = true;

if processStatic
    licorRetime = seconds(60);
    licorWindow = minutes(5);
end

if processDynamic
    daqRetime = seconds(30);
    daqWindow = minutes(30);
end

metadata = readtable("../DATA/FLUX/listing.xlsx");
licorDatasets = {};
daqDatasets = {};

fprintf("Importing Data:\n");

for currentDay = dateStart:dateEnd
    folderName = sprintf('%02d.%02d.%04d', currentDay.Month, currentDay.Day, currentDay.Year);
    dataDir = fullfile("../DATA/FLUX", folderName);
    
    if ~isfolder(dataDir)
        warning("Directory does not exist: %s", dataDir);
        continue;
    end
    
    files = dir(dataDir);
    
    for fileIdx = 1:length(files)
        file = files(fileIdx);
        
        if file.isdir
            continue;
        end
        
        rowIdx = strcmp(metadata{:,1}, file.name) & strcmp(metadata{:,2}, folderName);
        if ~any(rowIdx)
            warning("No metadata found for file: %s in folder: %s", file.name, folderName);
            continue;
        end
        fileMeta = metadata(rowIdx, :);
        filePath = fullfile(file.folder, file.name);
        
        if strcmpi(file.name, "licor.txt")
            if processStatic
                licor = importLicorData(filePath, fileMeta, cfg);
                if ~isempty(licor.data)
                    licorDatasets{end+1} = licor;
                end
            end
        elseif endsWith(file.name, ".txt", 'IgnoreCase', true)
            if processDynamic
                daq = importDAQData(filePath, fileMeta, startTime, endTime, cfg);
                if ~isempty(daq.data)
                    daqDatasets{end+1} = daq;
                end
            end
        end
    end
end

if ~processStatic
    licorDatasets = {};
end

if ~processDynamic
    daqDatasets = {};
end

processStatic = processStatic && ~isempty(licorDatasets);
processDynamic = processDynamic && ~isempty(daqDatasets);

licorCount = numel(licorDatasets);
daqCount = numel(daqDatasets);

fprintf("Imported %d LICOR datasets and %d DAQ datasets.\n", licorCount, daqCount);

if processStatic
    for i = 1:licorCount
        licor = licorDatasets{i};
        licor.data = retime(licor.data, "regular", 'mean', 'TimeStep', licorRetime);
        licor.data = varfun(@(x) smoothdata(x, 'movmean', licorWindow), licor.data);
        licor.data = unique(licor.data);
        licorDatasets{i} = licor;
    end
end

if processDynamic
    for i = 1:daqCount
        daq = daqDatasets{i};
        daq.data = retime(daq.data, "regular", 'mean', 'TimeStep', daqRetime);
        daq.data = smoothdata(daq.data, 'movmean', daqWindow);
        daq.data = unique(daq.data);
        daqDatasets{i} = daq;
    end
end

if processDynamic
    for i = 1:daqCount
        daq = daqDatasets{i};
        [daq.data, daq.meta] = autoCorrSensors(daq.data, daq.meta);
        daqDatasets{i} = daq;
    end
end

if processDynamic
    for i = 1:daqCount
        daq = daqDatasets{i};
        daq = applyDAQCalibration(daq, cfg);
        daqDatasets{i} = daq;
    end
end

if processStatic
    for i = 1:licorCount
        licor = licorDatasets{i};
        licor.data = calculateStaticFlux(licor.data, cfg);
        licorDatasets{i} = licor;
    end
end

if processDynamic
    for i = 1:daqCount
        daq = daqDatasets{i};
        daq.data = calculateDynamicFlux(daq.data, cfg);
        daqDatasets{i} = daq;
    end
end

expectedRange = cfg.mol_to_ppm([0.5e-6, 10e-6]);
expectedRange = cfg.ppm_to_mg(expectedRange);

licorBurn = [];
licorUnburn = [];
daqBurn = [];
daqUnburn = [];

if processStatic
    for i = 1:licorCount
        licor = licorDatasets{i};
        if strcmpi(licor.meta.DEPLOY_SITE, "BURN")
            licorBurn = [licorBurn; licor.data];
        else
            licorUnburn = [licorUnburn; licor.data];
        end
    end
    licorBurn = rmmissing(licorBurn);
    licorUnburn = rmmissing(licorUnburn);
end

if processDynamic
    for i = 1:daqCount
        daq = daqDatasets{i};
        if strcmpi(daq.meta.DEPLOY_SITE, "BURN")
            daqBurn = [daqBurn; daq.data(:, {"Q", "CA", "TA", "HA", "CB", "TB", "HB", "CA_CALIB", "CB_CALIB", "F_mg", "uF_mg"})];
        else
            daqUnburn = [daqUnburn; daq.data(:, {"Q", "CA", "TA", "HA", "CB", "TB", "HB", "CA_CALIB", "CB_CALIB", "F_mg", "uF_mg"})];
        end
    end
    daqBurn = rmmissing(daqBurn);
    daqUnburn = rmmissing(daqUnburn);
end

if processDynamic && ~isempty(daqBurn) && ~isempty(daqUnburn)
    daqBurnLumpDay = daqBurn.F_mg * 86400;
    daqUnburnLumpDay = daqUnburn.F_mg * 86400;
else
    daqBurnLumpDay = [];
    daqUnburnLumpDay = [];
end

if processStatic && ~isempty(licorBurn) && ~isempty(licorUnburn)
    licorBurnLumpDay = licorBurn.F_mg * 86400;
    licorUnburnLumpDay = licorUnburn.F_mg * 86400;
else
    licorBurnLumpDay = [];
    licorUnburnLumpDay = [];
end

fprintf("Reporting Data:\n");
fprintf("\tSITE\tSOURCE\tMEAN (mg/m^2/day)\tMEDIAN (mg/m^2/day)\n");

if ~isempty(daqBurnLumpDay)
    fprintf("\t%s\t%s\t%.4f\t\t%.4f\n", "Burn", "DAQ", mean(daqBurnLumpDay), median(daqBurnLumpDay));
end

if ~isempty(daqUnburnLumpDay)
    fprintf("\t%s\t%s\t%.4f\t\t%.4f\n", "Un-Burn", "DAQ", mean(daqUnburnLumpDay), median(daqUnburnLumpDay));
end

if ~isempty(licorBurnLumpDay)
    fprintf("\t%s\t%s\t%.4f\t\t%.4f\n", "Burn", "LICOR", mean(licorBurnLumpDay), median(licorBurnLumpDay));
end

if ~isempty(licorUnburnLumpDay)
    fprintf("\t%s\t%s\t%.4f\t\t%.4f\n", "Un-Burn", "LICOR", mean(licorUnburnLumpDay), median(licorUnburnLumpDay));
end

function licor = importLicorData(filePath, fileMeta, cfg)
    try
        licorData = IMPORTLICORFILE(filePath);
    catch ME
        warning("Failed to import LICOR file: %s. Error: %s", filePath, ME.message);
        licor = struct('data', table(), 'meta', table());
        return;
    end
    
    licorData.C_CALIB = licorData.C * 0.9883 + 24.6914;
    licor.meta = fileMeta;
    
    deploymentTimes = fileMeta{1, 15:22};
    dTimes = reshape(deploymentTimes, 2, [])';
    for dtIdx = 1:size(dTimes, 1)
        dRange = datetime(dTimes(dtIdx, :));
        if ~any(isnat(dRange))
            timeMask = (licorData.T > dRange(1)) & (licorData.T < dRange(2));
            tempData = licorData(timeMask, :);
            tempData = rmmissing(tempData);
            tempData.file = string(filePath);
            
            fprintf("\tFile Loaded: (LICOR)\t%s\n", filePath);
            
            switch dtIdx
                case {1, 2}
                    tempMeta.DEPLOY_SITE = "BURN";
                case {3, 4}
                    tempMeta.DEPLOY_SITE = "UNBURN";
            end
            licor.meta = tempMeta;
            
            if height(tempData) ~= 0
                licor.data = tempData;
            end
        end
    end
end

function daq = importDAQData(filePath, fileMeta, startTime, endTime, cfg)
    try
        daqData = IMPORTDAQFILE(filePath);
    catch ME
        warning("Failed to import DAQ file: %s. Error: %s", filePath, ME.message);
        daq = struct('data', table(), 'meta', table());
        return;
    end
    
    timeMask = (daqData.T > startTime) & (daqData.T < endTime);
    daqData = daqData(timeMask, :);
    daqData = rmmissing(daqData);
    
    daq.meta = fileMeta;
    daq.data = daqData;
    daq.file = string(filePath);
    
    fprintf("\tFile Loaded: (DAQ)\t%s\n", filePath);
end

function daq = applyDAQCalibration(daq, cfg)
    if strcmpi(daq.meta.PAD_A_ELT_SENSOR, "E") || strcmpi(daq.meta.PAD_B_ELT_SENSOR, "C")
        modelCA = load("../ELT CALIB/CALIB_W_LICOR/models/E-CL-linear-15-Jul-2024.mat");
        modelCB = load("../ELT CALIB/CALIB_W_LICOR/models/C-CL-linear-15-Jul-2024.mat");
        
        daq.data.CA_CALIB = predict(modelCA.model_2_lin, daq.data.CA);
        daq.data.CB_CALIB = predict(modelCB.model_1_lin, daq.data.CB);
    else
        daq.data.CA_CALIB = daq.data.CA;
        daq.data.CB_CALIB = daq.data.CB;
    end
end

function [data, meta] = autoCorrSensors(data, meta)
    caSTD = std(data.CA);
    cbSTD = std(data.CB);
    
    if caSTD > cbSTD
        [data.CA, data.CB] = deal(data.CB, data.CA);
        [data.TA, data.TB] = deal(data.TB, data.TA);
        [data.HA, data.HB] = deal(data.HB, data.HA);
    end
end

function data = calculateStaticFlux(data, cfg)
    dc = [0; diff(data.C_CALIB)];
    dt = [0; seconds(diff(data.T))];
    
    data.DCDT = dc ./ dt;
    data = rmmissing(data);
    
    data.DCDT_mg = cfg.ppm_to_mg(data.DCDT);
    data.C_mg = cfg.ppm_to_mg(data.C_CALIB);
    data.F_mg = (data.DCDT_mg .* cfg.V) ./ cfg.As;
end

function data = calculateDynamicFlux(data, cfg)
    Q = cfg.lpm_to_cms(mean(data.Q / 1000));
    
    if Q == 0
        Q = cfg.lpm_to_cms(2.0204125);
    end
    
    AMBIENT = cfg.ppm_to_mg(data.CA_CALIB);
    CHAMBER = cfg.ppm_to_mg(data.CB_CALIB);
    CHAMBER_FLOOR = CHAMBER - AMBIENT(1);
    
    data.F_mg = (Q .* (CHAMBER - AMBIENT)) ./ cfg.As;
    
    dC = cfg.ppm_to_mg(data.CB_CALIB - data.CA_CALIB);
    uQ = cfg.lpm_to_cms(cfg.uQ);
    uC = cfg.ppm_to_mg(cfg.modelCB.model_1_lin.RMSE);
    
    data.uF_mg = sqrt((dC .* uQ ./ cfg.As).^2 + ((2 .* Q .* uC) ./ cfg.As).^2);
end
