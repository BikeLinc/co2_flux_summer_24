%% calibration.m
% This script generates calibration models for all sensors listed in the CSV.

clc;
clear;
close all;

if ~isfolder('figs_training')
    mkdir('figs_training');
end

addpath(genpath('../../UTILS'));

%% Importing Data
disp("Importing Data");
csvFilePath = 'overlapping_periods.csv';
fileMapping = readtable(csvFilePath);
disp(['Number of rows found in CSV file: ', num2str(height(fileMapping))]);

selection = questdlg('Select data type to include:', 'Data Selection', ...
                     'AMB', 'CL', 'Both', 'Both');
if isempty(selection)
    error('No selection made. Script terminated.');
end

selectedFiles = filterFileMapping(fileMapping, selection);
disp(['Number of files selected: ', num2str(height(selectedFiles))]);

%% Initialize
[daqData, refData, settings] = initialize();

save('cache/temp_daq_ref.mat', 'daqData', 'refData', 'settings');
clear daqData refData;

%% Importing Data
load('cache/temp_daq_ref.mat', 'daqData', 'refData', 'settings');
[daqData, refData] = importSelectedFiles(selectedFiles, daqData, refData, settings);

save('cache/temp_imported_data.mat', 'daqData', 'refData');
clear daqData refData;

%% Collect each sensor’s data
load('cache/temp_imported_data.mat', 'daqData', 'refData');
allSensors = batchSensorData(selectedFiles, daqData, refData, settings);

save('cache/temp_all_sensors.mat', 'allSensors');
clear allSensors;

%% Clean Sensor Data
load('cache/temp_all_sensors.mat', 'allSensors');
allSensors = cleanSensorData(allSensors, settings.errorValues);

save('cache/temp_cleaned_sensors.mat', 'allSensors');
clear allSensors;

%% Plot Raw Data
load('cache/temp_cleaned_sensors.mat', 'allSensors');
plotRawData(allSensors, refData, "figs_training/raw_data_overlap.fig");

save('cache/temp_plot_data.mat', 'allSensors');
clear allSensors;

%% Further Pre-Processing and Synchronization
load('cache/temp_plot_data.mat', 'allSensors');
allSensors = furtherPreprocessAndSync(allSensors, refData, settings);

save('cache/temp_synced_sensors.mat', 'allSensors');
clear allSensors;

%% Parameter Definitions
[params, predictors, targetVariable] = defineParameters();

%% Data Replication
load('cache/temp_synced_sensors.mat', 'allSensors');
allSensors = replicateData(allSensors, predictors, params);

save('cache/temp_replicated_data.mat', 'allSensors');
clear allSensors;

%% Partition Data into Train/Validation/Evaluation Sets
load('cache/temp_replicated_data.mat', 'allSensors');
dataPartitions = partitionData(allSensors, targetVariable, params);

save('cache/temp_data_partitions.mat', 'dataPartitions');
clear dataPartitions;

%% Train Models (Linear & Neural Network)
load('cache/temp_data_partitions.mat', 'dataPartitions');
[models, modelComparisons] = trainModels(allSensors, dataPartitions, predictors, targetVariable, params);

save('cache/temp_models.mat', 'models', 'modelComparisons');
clear models modelComparisons;

%% Create and Display Overall Model Comparison Table
load('cache/temp_models.mat', 'modelComparisons');
createOverallModelComparison(modelComparisons, "model_comparison_all_sensors.csv");

clear modelComparisons;

%% Export Models
load('cache/temp_models.mat', 'models');
exportModels(models, "models", datestr(datetime("today"), 'yyyymmdd'));

clear models;
disp("All requested sensor calibrations complete.");

%% Evaluate Other Models

disp("Evaluating previously saved models on a new dataset...");

evaluateSavedModels(dataPartitions.E.evaluationData, ["X_C","X_T","X_H"], "Y_C");

%% Function Definitions

function selectedFiles = filterFileMapping(fileMapping, selection)
    switch selection
        case 'AMB'
            selectedFiles = fileMapping(strcmpi(fileMapping.Type, 'AMB'), :);
        case 'CL'
            selectedFiles = fileMapping(strcmpi(fileMapping.Type, 'CL'), :);
        case 'Both'
            selectedFiles = fileMapping;
        otherwise
            error('Invalid selection. Script terminated.');
    end
end

function [daqData, refData, settings] = initialize()
    daqData = [];
    refData = [];

    settings.smoothDuration      = minutes(15);
    settings.retimeDuration      = minutes(1);
    settings.outlierPercentiles  = [2, 98]; 
    settings.outlierRemoval      = false;
    settings.errorValues         = [-999, 500, 2815, 64537, 231753, 65535, 2500, 2559];
    settings.referenceMin        = 390;
    settings.referenceMax        = 1200;
end

function [daqData, refData] = importSelectedFiles(selectedFiles, daqData, refData, settings)
    disp("Pre-Processing Data");
    for i = 1:height(selectedFiles)
        daqFilePath     = "../../DATA/CALIB/" + selectedFiles.DAQFile{i};
        licorFilePath   = "../../DATA/CALIB/" + selectedFiles.LicorFile{i};
        picarroFilePath = "../../DATA/CALIB/" + selectedFiles.PICARROFile{i};
    
        overlapStart = datetime(selectedFiles.OverlapStart(i), 'InputFormat','MM/dd/yyyy HH:mm:ss');
        overlapEnd   = datetime(selectedFiles.OverlapEnd(i), 'InputFormat','MM/dd/yyyy HH:mm:ss');
    
        daqDataTmp = IMPORTDAQFILE(daqFilePath);
        daqDataTmp = daqDataTmp(:, ["CA", "TA", "HA", "CB", "TB", "HB"]);
    
        if isempty(daqDataTmp) || ~istimetable(daqDataTmp)
            disp(['DAQ data is empty or invalid for file: ', daqFilePath, '. Skipping...']);
            continue;
        end
    
        try
            if contains(licorFilePath, 'licor')
                referenceDataTmp = IMPORTLICORFILE(licorFilePath);
                referenceDataTmp = referenceDataTmp(referenceDataTmp.C >= settings.referenceMin & referenceDataTmp.C <= settings.referenceMax, :);
                referenceDataTmp = referenceDataTmp(:, "C");
            else
                referenceDataTmp = readtable(picarroFilePath, 'FileType', 'text');
                
                if ismember("CO2_sync", referenceDataTmp.Properties.VariableNames)
                    referenceDataTmp = referenceDataTmp(referenceDataTmp.CO2_sync >= settings.referenceMin & referenceDataTmp.CO2_sync <= settings.referenceMax, :);
                    referenceDataTmp = referenceDataTmp(:, ["T", "CO2_sync"]);
                else
                    referenceDataTmp = referenceDataTmp(referenceDataTmp.CO2 >= settings.referenceMin & referenceDataTmp.CO2 <= settings.referenceMax, :);
                    referenceDataTmp = referenceDataTmp(:, ["T", "CO2"]);
                end
                referenceDataTmp.Properties.VariableNames = ["T", "C"];
                referenceDataTmp.T = referenceDataTmp.T - hours(7);
                referenceDataTmp = table2timetable(referenceDataTmp);
            end
        catch
            if licorFilePath ~= "../../DATA/CALIB/"
                disp(['Error reading reference data for file: ', char(licorFilePath), '. Skipping...']);
            else
                disp(['Error reading reference data for file: ', char(picarroFilePath), '. Skipping...']);
            end
            continue;
        end
    
        referenceDataCurrent = referenceDataTmp;
    
        if ~isempty(daqDataTmp) 
            daqDataTmp = removeErrorValues(daqDataTmp, selectedFiles, i, settings.errorValues);
    
            daqDataTmp = removeDuplicatesAndMissing(daqDataTmp);
            referenceDataCurrent = removeDuplicatesAndMissing(referenceDataCurrent);

            daqDataTmp = sortrows(daqDataTmp);
            referenceDataCurrent = sortrows(referenceDataCurrent);

            daqDataTmp = retime(daqDataTmp, 'regular','mean','TimeStep', settings.retimeDuration);
            referenceDataCurrent = retime(referenceDataCurrent, 'regular','mean','TimeStep', settings.retimeDuration);

            [referenceDataCurrent, daqDataTmp] = alignDataWithCrossCorr(referenceDataCurrent, daqDataTmp);
    
            if settings.outlierRemoval
                daqDataTmp = rmoutliers(daqDataTmp, 'percentile', settings.outlierPercentiles);
                referenceDataCurrent = rmoutliers(referenceDataCurrent, 'percentile', settings.outlierPercentiles);
            end
    
            daqDataTmp = removeMissingRows(daqDataTmp, 'T');
            referenceDataCurrent = removeMissingRows(referenceDataCurrent, 'T');            
    
            daqDataTmp = smoothTimetableData(daqDataTmp, settings.smoothDuration);
            referenceDataCurrent = smoothTimetableData(referenceDataCurrent, settings.smoothDuration);
    
            disp("+ " + max(height(daqDataTmp), height(referenceDataCurrent)) + " lines");
    
            daqData = [daqData; daqDataTmp];
            refData = [refData; referenceDataCurrent];
        end
    end

    daqData = removeDuplicatesAndMissing(daqData);
    refData = removeDuplicatesAndMissing(refData);
    
    daqMissingRowsFinal = any(ismissing(daqData), 2);
    daqData(daqMissingRowsFinal, :) = [];
    
    disp("NaN Rows Removed: " + sum(daqMissingRowsFinal));
end

function daqDataTmp = removeErrorValues(daqDataTmp, selectedFiles, i, errorValues)
    if ismember('CA', daqDataTmp.Properties.VariableNames) && selectedFiles.SensorA{i} ~= "X"
        daqDataTmp(ismember(daqDataTmp.CA, errorValues), :) = [];
        daqDataTmp = rmmissing(daqDataTmp, "DataVariables", ["CA","TA","HA"]);
    end
    if ismember('CB', daqDataTmp.Properties.VariableNames) && selectedFiles.SensorB{i} ~= "X"
        daqDataTmp(ismember(daqDataTmp.CB, errorValues), :) = [];
        daqDataTmp = rmmissing(daqDataTmp, "DataVariables", ["CB","TB","HB"]);
    end
end

function timetableOut = removeDuplicatesAndMissing(timetableIn)
    [~, uniqueIdx] = unique(timetableIn.T);
    timetableOut = timetableIn(uniqueIdx, :);
    timetableOut = rmmissing(timetableOut);
end

function timetableOut = removeMissingRows(timetableIn, timeVar)
    missingRows = any(ismissing(timetableIn), 2) | isnat(timetableIn.(timeVar));
    timetableOut = timetableIn;
    timetableOut(missingRows, :) = [];
end

function timetableOut = smoothTimetableData(timetableIn, smoothDuration)
    smoothDurationCalc = round(height(timetableIn)*0.01);
    if smoothDurationCalc > height(timetableIn)
        smoothDurationCalc = height(timetableIn);
    end
    if smoothDurationCalc < 1
        smoothDurationCalc = 1;
    end
    timetableOut = timetableIn;
    vars = timetableIn.Properties.VariableNames;
    for v = 1:length(vars)
        if isnumeric(timetableIn.(vars{v}))
            timetableOut.(vars{v}) = smoothdata(timetableIn.(vars{v}), 'movmean', smoothDurationCalc);
        end
    end
end

function allSensors = batchSensorData(selectedFiles, daqData, refData, settings)
    disp("--> Separating Sensor Data for ALL sensors");
    allSensors = struct();
    
    for i = 1:height(selectedFiles)
        currentSensorA = selectedFiles.SensorA{i};
        currentSensorB = selectedFiles.SensorB{i};
        startTime      = datetime(selectedFiles.OverlapStart(i), 'InputFormat','MM/dd/yyyy HH:mm:ss');
        endTime        = datetime(selectedFiles.OverlapEnd(i), 'InputFormat','MM/dd/yyyy HH:mm:ss');
    
        idxDaq = daqData.T >= startTime & daqData.T <= endTime;
        currDaq = daqData(idxDaq, :);
        
        if isempty(currDaq)
            disp(['No DAQ data found for overlap period: ', datestr(startTime), ' to ', datestr(endTime), '. Skipping...']);
            continue;
        end
    
        A_data = currDaq(:, ["CA", "TA", "HA"]);
        B_data = currDaq(:, ["CB", "TB", "HB"]);
        
        A_data.Properties.VariableNames = ["CA","TA","HA"];
        B_data.Properties.VariableNames = ["CA","TA","HA"];
    
        if ~isfield(allSensors, currentSensorA)
            allSensors.(currentSensorA) = timetable();
        end
        allSensors.(currentSensorA) = [allSensors.(currentSensorA); A_data];
    
        if ~isfield(allSensors, currentSensorB)
            allSensors.(currentSensorB) = timetable();
        end
        allSensors.(currentSensorB) = [allSensors.(currentSensorB); B_data];
    end
end

function allSensors = cleanSensorData(allSensors, errorValues)
    disp("Cleaning Sensor Data");
    
    fieldsList = fieldnames(allSensors);
    for f = 1:numel(fieldsList)
        sensorID  = fieldsList{f};
        sensorTable = allSensors.(sensorID);
        
        [~, uniqueIdx] = unique(sensorTable.T, 'stable');
        sensorTable = sensorTable(uniqueIdx, :);
        sensorTable = rmmissing(sensorTable);
        
        sensorTable(ismember(sensorTable.CA, errorValues), :) = [];
        
        if isempty(sensorTable)
            disp(['No data available for sensor ', sensorID]);
            continue;
        end
        
        allSensors.(sensorID) = sensorTable;
    end
end

function plotRawData(allSensors, refData, savePath)
    disp("Plotting Raw Data for Reference (All Sensors)");
    
    fieldsList = fieldnames(allSensors);
    fig = figure('Name', 'Raw Data Overlap - All Sensors', 'Position', [100 100 1200 800]); hold on; grid on;
    plot(refData.T, refData.C, 'Marker','o','LineStyle','none','DisplayName','LICOR');
    
    for f = 1:numel(fieldsList)
        sensorID  = fieldsList{f};
        sensorTable = allSensors.(sensorID);
        if ~isempty(sensorTable)
            plot(sensorTable.T, sensorTable.CA, 'Marker','.','LineStyle','none','DisplayName',['Sensor ' sensorID]);
        end
    end
    
    xlabel("Time"); ylabel("CO_2 [ppm]");
    title("Raw Data Overlap - All Sensors");
    legend('Location','eastoutside');
    
    saveas(fig, savePath);
    close all;
end

function allSensors = furtherPreprocessAndSync(allSensors, refData, settings)
    disp("Further Pre-Processing and Synchronization");
    fieldsList = fieldnames(allSensors);
    
    for f = 1:numel(fieldsList)
        sensorID  = fieldsList{f};
        sensorTable = allSensors.(sensorID);
    
        if ~isempty(sensorTable)
            if settings.outlierRemoval
                sensorTable = rmoutliers(sensorTable, 'percentile', settings.outlierPercentiles);
            end
            
            sensorTable = retime(sensorTable, "regular", 'mean','TimeStep', settings.retimeDuration);
        
            sensorTable = synchronize(sensorTable, refData);
            sensorTable = rmmissing(sensorTable);
        
            sensorTable = renamevars(sensorTable, [1,2,3,4], ["X_C","X_T","X_H","Y_C"]);
        
            allSensors.(sensorID) = sensorTable;
        end
    end
end

function [params, predictors, targetVariable] = defineParameters()
    disp("Parameter Definitions");
    
    params.numBins = 25;
    params.targetBinCount = 150;
    params.trainFraction = 0.50;
    params.validationFraction = 0.1;
    params.evaluationFraction = 0.49;
    params.replicateData = true;
    params.useRandomSplit = true;
    
    predictors = ["X_C", "X_T", "X_H"];  
    targetVariable = "Y_C";
end

function allSensors = replicateData(allSensors, predictors, params)
    disp("Data Replication (if enabled)");
    
    fieldsList = fieldnames(allSensors);
    for f = 1:numel(fieldsList)
        sensorID = fieldsList{f};
        if sensorID == "X", continue; end
    
        sensorTable = timetable2table(allSensors.(sensorID));
        if isempty(sensorTable), continue; end
    
        if params.replicateData
            disp("Replicating data for sensor " + sensorID + " ...");
    
            xData = sensorTable{:, predictors(1)};
            minVal = min(xData);
            maxVal = max(xData);
    
            if minVal == maxVal
                disp("No variation in " + string(predictors(1)) + ", skipping replication.");
                continue;
            end
    
            binEdges  = linspace(minVal, maxVal, params.numBins+1);
            binIdx    = discretize(xData, binEdges);
    
            replicatedTable = table();
            for b = 1:params.numBins
                binMask = (binIdx == b);
                binData = sensorTable(binMask, :);
                if isempty(binData), continue; end
    
                binCount = height(binData);
                if binCount < params.targetBinCount
                    repsNeeded = params.targetBinCount - binCount;
                    repIdx     = randi([1, binCount], repsNeeded, 1);
                    binData    = [binData; binData(repIdx, :)];
                end
                replicatedTable = [replicatedTable; binData];
            end
            sensorTable = replicatedTable;
        end
    
        if ~isempty(sensorTable)
            allSensors.(sensorID) = table2timetable(sensorTable, 'RowTimes', sensorTable.T);
        end
    end
end

function dataPartitions = partitionData(allSensors, targetVariable, params)
    disp("Partitioning Data into Train/Validation/Evaluation Sets");
    
    dataPartitions = struct();
    fieldsList = fieldnames(allSensors);
    
    for f = 1:numel(fieldsList)
        sensorID = fieldsList{f};
        if sensorID == "X", continue; end
    
        sensorTable = timetable2table(allSensors.(sensorID));
        if height(sensorTable) < 5
            disp("Too few data points for sensor " + sensorID + "; skipping partition.");
            continue;
        end
    
        targetData = sensorTable.(targetVariable);
        targetMin  = min(targetData);
        targetMax  = max(targetData);
    
        if targetMin == targetMax
            disp("No range in " + targetVariable + " for sensor " + sensorID + "; skipping partition.");
            continue;
        end
    
        binEdges = linspace(targetMin, targetMax, params.numBins+1);
    
        trainData = table();
        validationData = table();
        evaluationData = table();
    
        for b = 1:params.numBins
            inBin = (targetData >= binEdges(b) & targetData < binEdges(b+1));
            binData = sensorTable(inBin, :);
            if isempty(binData), continue; end
    
            if height(binData) < params.targetBinCount, continue; end
    
            binData  = sortrows(binData, 'Time', 'ascend');
    
            nBin     = height(binData);
            trainCut = floor(params.trainFraction * nBin);
            valCut   = floor((params.trainFraction + params.validationFraction) * nBin);
    
            if params.useRandomSplit
                binData = binData(randperm(nBin), :);
            end
    
            binTrainData = binData(1 : trainCut, :);
            binValidationData = binData(trainCut+1 : valCut, :);
            binEvaluationData = binData(valCut+1 : end, :);
    
            trainData = [trainData; binTrainData];
            validationData = [validationData; binValidationData];
            evaluationData = [evaluationData; binEvaluationData];
        end
    
        disp("Sensor: " + sensorID);
        disp("   Training Data:   " + num2str(height(trainData)) + " rows.");
        disp("   Validation Data: " + num2str(height(validationData))   + " rows.");
        disp("   Evaluation Data: " + num2str(height(evaluationData))  + " rows.");
    
        dataPartitions.(sensorID).trainData = trainData;
        dataPartitions.(sensorID).validationData = validationData;
        dataPartitions.(sensorID).evaluationData  = evaluationData;
    
        plotPartitions(sensorID, trainData, validationData, evaluationData);
    end
    
    function plotPartitions(sensorID, trainData, validationData, evaluationData)
        fig = figure('Name', "Partition: " + sensorID, 'Position', [200 200 1200 500]);
        tiledlayout(3,1, 'Padding', 'compact');
    
        nexttile; hold on; grid on;
        plot(trainData.Time, trainData.Y_C, 'b.');
        title("Training Data"); xlabel("Time"); ylabel("Y_C");
    
        nexttile; hold on; grid on;
        plot(validationData.Time, validationData.Y_C, 'g.');
        title("Validation Data"); xlabel("Time"); ylabel("Y_C");
    
        nexttile; hold on; grid on;
        plot(evaluationData.Time, evaluationData.Y_C, 'r.');
        title("Evaluation Data"); xlabel("Time"); ylabel("Y_C");
    
        saveas(fig, "figs_training/PARTITION_DATA_"+sensorID+".fig")
        close all;
    end
end

function [models, modelComparisons] = trainModels(allSensors, dataPartitions, predictors, targetVariable, params)
    disp("Training Models (Linear & Neural Network)");
    waitbr = waitbar(0, "Training Models...");

    allCombinations = {{'X_C','X_T','X_H'}};

    referenceRanges = { [390, 450], [450, 1200], [390, 1200] };
    rangeLabels     = { '390_450',  '450_1200',  '390_1200'  };

    models = struct();  
    modelComparisons = struct();
    fieldsList = fieldnames(allSensors);

    for f = 1:numel(fieldsList)
        sensorID = fieldsList{f};
        if sensorID == "X", continue; end

        if ~isfield(dataPartitions, sensorID) || ~isfield(allSensors, sensorID)
            disp("No data available for sensor " + sensorID + "; skipping training.");
            continue;
        end

        allData = timetable2table(allSensors.(sensorID));

        modelComparison = table('Size', [0, 6], ...
            'VariableTypes', {'string','string','string','string','double','double'}, ...
            'VariableNames', {'SensorID','Range','ModelType','Predictors','RMSE','R2'});

        for rr = 1:numel(referenceRanges)
            currentRange = referenceRanges{rr};
            rangeLabel   = rangeLabels{rr};

            allDataSubset = allData(allData.(targetVariable) >= currentRange(1) & ...
                                    allData.(targetVariable) <= currentRange(2), :);

            trainData = filterRange(dataPartitions.(sensorID).trainData, targetVariable, currentRange);
            validationData = filterRange(dataPartitions.(sensorID).validationData, targetVariable, currentRange);
            evaluationData = filterRange(dataPartitions.(sensorID).evaluationData, targetVariable, currentRange);

            if isempty(trainData) || isempty(validationData) || isempty(evaluationData) || isempty(allDataSubset)
                disp("Not enough data in [" + currentRange(1) + "," + currentRange(2) + ...
                     "] for sensor " + sensorID + "; skipping...");
                continue;
            end

            fig = figure();
            hold on; grid on;
            reference = [evaluationData; validationData; trainData];
            if rangeLabel == "390_450"
                save("../../PARAM_SPACE/amb_calibration_param.mat", 'reference')
            elseif rangeLabel == "450_1200"
                save("../../PARAM_SPACE/cl_calibration_param.mat", 'reference')
            else
                save("../../PARAM_SPACE/comb_calibration_param.mat", 'reference')
            end
            plot(reference.T, reference.Y_C, 'k.', 'DisplayName', 'Reference');
            plot(reference.T, reference.X_C, 'b.', 'DisplayName', 'LCNDIR');
            legend();
            xlabel("Time")
            ylabel("CO_2 (ppm)")
            saveas(fig, "figs_training/SENSOR_" + sensorID + "_DATA_" + rangeLabel + ".fig");
            close all;

            for c = 1:length(allCombinations)
                waitbar(c/length(allCombinations), waitbr, "Training Models...");
                currentPredictors = allCombinations{c};

                predictorStr_display = strjoin(currentPredictors, ',');
                predictorStr_struct  = strjoin(currentPredictors, '_');

                try
                    linModel = fitlm(allDataSubset{:, currentPredictors}, ...
                                     allDataSubset.(targetVariable), 'interactions');
                catch ME
                    disp("Error training Linear Model for Sensor " + sensorID + ...
                         " with predictors [" + predictorStr_display + "]: " + ME.message);
                    continue;
                end

                try
                    linPred_all = predict(linModel, allDataSubset{:, currentPredictors});
                    y_all       = allDataSubset.(targetVariable);
                    linRMSE_all = sqrt(mean((linPred_all - y_all).^2, 'omitnan'));
                    ssRes_lin   = sum((y_all - linPred_all).^2, 'omitnan');
                    ssTot_lin   = sum((y_all - mean(y_all)).^2, 'omitnan');
                    linR2_all   = 1 - ssRes_lin/ssTot_lin;
                catch ME
                    disp("Error computing metrics for Linear Model for Sensor " + sensorID + ...
                         " with predictors [" + predictorStr_display + "]: " + ME.message);
                    linRMSE_all = NaN; linR2_all = NaN;
                end

                if isempty(trainData) || isempty(validationData) || isempty(evaluationData)
                    disp("Missing train/validation/evaluation data for sensor " + sensorID + "; skipping Neural Net.");
                    continue;
                end

                maxTries     = 10;
                targetRMSE   = 15;
                bestNet      = [];
                bestValRMSE  = inf;

                for attempt = 1:maxTries
                    try
                        tempNet = feedforwardnet([16,16]);
                        tempNet.trainParam.epochs = 1000;
                        tempNet.trainParam.showWindow = false;
                        tempNet = train(tempNet, trainData{:, currentPredictors}', ...
                                                trainData.(targetVariable)');
                    catch ME
                        disp("Error training NN for Sensor " + sensorID + " [" + predictorStr_display + ...
                             "]: " + ME.message);
                        break;
                    end

                    try
                        valNetPred = tempNet(validationData{:, currentPredictors}')';
                        valY       = validationData.(targetVariable);
                        valNetRMSE = sqrt(mean((valNetPred - valY).^2, 'omitnan'));
                    catch ME
                        disp("Error computing validation RMSE for Sensor " + sensorID + ...
                             " [" + predictorStr_display + "]: " + ME.message);
                        continue;
                    end

                    if valNetRMSE < bestValRMSE
                        bestValRMSE = valNetRMSE;
                        bestNet     = tempNet;
                    end

                    if valNetRMSE < targetRMSE
                        disp("Sensor " + sensorID + ": validation RMSE < " + targetRMSE + ...
                             " with [" + predictorStr_display + "] on attempt " + attempt);
                        break;
                    else
                        disp("Sensor " + sensorID + ": attempt " + attempt + ...
                             " -> val RMSE = " + valNetRMSE + " (target < " + targetRMSE + ")");
                    end
                end

                if bestValRMSE < targetRMSE
                    disp("Sensor " + sensorID + ": final NN meets RMSE < " + targetRMSE + ...
                         " with [" + predictorStr_display + "]");
                else
                    disp("Sensor " + sensorID + ": NN did not reach RMSE < " + targetRMSE + ...
                         " after " + maxTries + " attempts. Using best so far...");
                end

                rangeKey = [predictorStr_struct '_' rangeLabel];
                models.(sensorID).(rangeKey).linearModel = linModel;
                models.(sensorID).(rangeKey).neuralNet   = bestNet;

                try
                    valX = validationData{:, currentPredictors};
                    valY = validationData.(targetVariable);
                    valLinPred = predict(linModel, valX);
                    valLinRMSE = sqrt(mean((valLinPred - valY).^2, 'omitnan'));
                    ssRes      = sum((valY - valLinPred).^2, 'omitnan');
                    ssTot      = sum((valY - mean(valY)).^2, 'omitnan');
                    valLinR2   = 1 - ssRes/ssTot;

                    valNetPred = bestNet(valX')';
                    valNetRMSE = sqrt(mean((valNetPred - valY).^2, 'omitnan'));
                    ssResN     = sum((valY - valNetPred).^2, 'omitnan');
                    ssTotN     = sum((valY - mean(valY)).^2, 'omitnan');
                    valNetR2   = 1 - ssResN/ssTotN;
                catch ME
                    disp("Error computing validation metrics for Sensor " + sensorID + ...
                         " [" + predictorStr_display + "]: " + ME.message);
                    valLinRMSE = NaN; valLinR2 = NaN;
                    valNetRMSE = NaN; valNetR2 = NaN;
                end

                try
                    evalX = evaluationData{:, currentPredictors};
                    evalY = evaluationData.(targetVariable);

                    evalLinPred  = predict(linModel, evalX);
                    evalLinRMSE  = sqrt(mean((evalLinPred - evalY).^2, 'omitnan'));
                    ssRes_eval_l = sum((evalY - evalLinPred).^2, 'omitnan');
                    ssTot_eval_l = sum((evalY - mean(evalY)).^2, 'omitnan');
                    evalLinR2    = 1 - ssRes_eval_l/ssTot_eval_l;

                    if isempty(bestNet)
                        evalNetPred = NaN(size(evalY));
                        evalNetRMSE = NaN; 
                        evalNetR2   = NaN;
                    else
                        evalNetPred = bestNet(evalX')';
                        evalNetRMSE = sqrt(mean((evalNetPred - evalY).^2, 'omitnan'));
                        ssRes_eval_n = sum((evalY - evalNetPred).^2, 'omitnan');
                        ssTot_eval_n = sum((evalY - mean(evalY)).^2, 'omitnan');
                        evalNetR2    = 1 - ssRes_eval_n/ssTot_eval_n;
                    end
                catch ME
                    disp("Error computing evaluation metrics: " + ME.message);
                    evalLinRMSE = NaN; evalLinR2 = NaN; 
                    evalNetRMSE = NaN; evalNetR2 = NaN;
                end

                disp("Final Evaluation [" + rangeLabel + "] (Sensor " + sensorID + ...
                     ", Predictors [" + predictorStr_display + "])");
                disp("Linear  Model - RMSE: " + evalLinRMSE + ", R^2: " + evalLinR2);
                disp("Neural Network Model - RMSE: " + evalNetRMSE + ", R^2: " + evalNetR2);

                plotResiduals(sensorID + "_" + rangeLabel, predictorStr_display, evalY, ...
                              evalLinPred, evalNetPred, targetVariable, ...
                              evalLinRMSE, evalLinR2, evalNetRMSE, evalNetR2);

                plotTimeSeries(sensorID + "_" + rangeLabel, predictorStr_display, evaluationData.Time, ...
                               evalY, evalLinPred, evalNetPred, ...
                               evalLinRMSE, evalLinR2, evalNetRMSE, evalNetR2);

                modelComparison = [modelComparison; 
                    {sensorID, rangeLabel, 'Linear', predictorStr_display, evalLinRMSE, evalLinR2}; 
                    {sensorID, rangeLabel, 'Neural Network', predictorStr_display, evalNetRMSE, evalNetR2}];
            end
        end

        modelComparisons.(sensorID) = modelComparison;
        disp("Model Comparison for Sensor " + sensorID);
        disp(modelComparison);

        save("models_created_checkpoint.mat");
    end

    close(waitbr);

    function subData = filterRange(tbl, tgtVar, rangeVec)
        if isempty(tbl)
            subData = tbl;
            return;
        end
        subData = tbl(tbl.(tgtVar) >= rangeVec(1) & tbl.(tgtVar) <= rangeVec(2), :);
    end

    function plotResiduals(sensorTag, predictorStr_display, trueY, linPred, evalNetPred, ...
                           targetVariable, linRMSE, linR2, netRMSE, netR2)
        fig = figure('Name', "Residual Plots - Sensor " + sensorTag + ...
                           " - Predictors [" + predictorStr_display + "]", ...
                     'Position', [100 100 1200 800]);
        tiledlayout(2,2, 'Padding', 'compact');

        nexttile; hold on; grid on;
        scatter(trueY, linPred - trueY, 'b.');
        yline(0, 'r--','HandleVisibility','off');
        xlabel("True " + targetVariable); ylabel("Residual (Linear)");
        title("Eval Residuals - Linear Model");
        legend(sprintf('Residuals\nRMSE: %.2f, R²: %.2f', linRMSE, linR2), 'Location', 'best');

        nexttile; hold on; grid on;
        scatter(trueY, evalNetPred - trueY, 'g.');
        yline(0, 'r--','HandleVisibility','off');
        xlabel("True " + targetVariable); ylabel("Residual (Neural Network)");
        title("Eval Residuals - Neural Network");
        legend(sprintf('Residuals\nRMSE: %.2f, R²: %.2f', netRMSE, netR2), 'Location','best');

        nexttile; hold on; grid on;
        scatter(trueY, linPred, 'b.');
        plot(trueY, trueY, 'r--','HandleVisibility','off');
        xlabel("True " + targetVariable);
        ylabel("Prediction");
        title("Eval 1:1 - Linear Model");
        legend(sprintf('Predictions\nRMSE: %.2f, R²: %.2f', linRMSE, linR2), 'Location','best');

        nexttile; hold on; grid on;
        scatter(trueY, evalNetPred, 'g.');
        plot(trueY, trueY, 'r--','HandleVisibility','off');
        xlabel("True " + targetVariable); ylabel("Prediction");
        title("Time Series - Neural Network Predictions");
        legend(sprintf('Predictions\nRMSE: %.2f, R²: %.2f', netRMSE, netR2), 'Location','best');

        saveas(fig, "figs_training/SENSOR_" + sensorTag + "_PREDICTORS_[" + predictorStr_display + "].fig");
        close all;
    end

    function plotTimeSeries(sensorTag, predictorStr_display, timeVec, trueY, linPred, evalNetPred, ...
                            linRMSE, linR2, netRMSE, netR2)
        if isrow(timeVec),       timeVec = timeVec';  end
        if isrow(trueY),         trueY = trueY';      end
        if isrow(linPred),       linPred = linPred';  end
        if isrow(evalNetPred),   evalNetPred = evalNetPred'; end

        if ~(length(timeVec) == length(trueY) && ...
             length(linPred) == length(evalNetPred) && ...
             length(linPred) == length(trueY))
            error('Time vector and predictions must be same length.');
        end

        disp("Generating Time Series Plots for Sensor " + sensorTag + ...
             " with Predictors [" + predictorStr_display + "]");

        fig = figure('Name', "Time Series - Sensor " + sensorTag + ...
                           " - Predictors [" + predictorStr_display + "]", ...
                     'Position', [100 100 1600 900]);
        tiledlayout(2,1, 'Padding', 'compact');

        nexttile; hold on; grid on;
        plot(timeVec, trueY, 'k.', 'DisplayName', 'True Y_C');
        plot(timeVec, linPred, 'b.', 'DisplayName', ...
             sprintf('Linear Model\nRMSE: %.2f, R²: %.2f', linRMSE, linR2));
        xlabel("Time"); ylabel("Y_C");
        title("Time Series - Linear Model Predictions");
        legend('Location','best');

        nexttile; hold on; grid on;
        plot(timeVec, trueY, 'k.', 'DisplayName', 'True Y_C');
        plot(timeVec, evalNetPred, 'r.', 'DisplayName', ...
             sprintf('Neural Network\nRMSE: %.2f, R²: %.2f', netRMSE, netR2));
        xlabel("Time"); ylabel("Y_C");
        title("Time Series - Neural Network Predictions");
        legend('Location','best');

        saveas(fig, "figs_training/SENSOR_" + sensorTag + "_TIMESERIES_[" + predictorStr_display + "].fig");
        close all;
    end
end

function createOverallModelComparison(modelComparisons, savePath)
    disp("Creating Overall Model Comparison Table");
    
    overallModelComparison = table();
    fieldsList = fieldnames(modelComparisons);
    for f = 1:numel(fieldsList)
        sensorID = fieldsList{f};
        if ~isfield(modelComparisons, sensorID), continue; end
        sensorModels = modelComparisons.(sensorID);
        overallModelComparison = [overallModelComparison; sensorModels];
    end
    
    disp("Overall Model Comparison");
    disp(sortrows(overallModelComparison, 'R2', 'descend'));
    
    writetable(overallModelComparison, savePath);
end

function exportModels(models, folderName, modelSaveDate)
    disp("Exporting Final Models");
    
    if ~isfolder(folderName)
        mkdir(folderName);
    end
    fieldsList = fieldnames(models);
    for f = 1:numel(fieldsList)
        sensorID = fieldsList{f};
        if ~isfield(models, sensorID), continue; end
        predictorSets = fieldnames(models.(sensorID));
        for p = 1:numel(predictorSets)
            predictorSet = predictorSets{p};
            linearModel = models.(sensorID).(predictorSet).linearModel;
            neuralNetModel = models.(sensorID).(predictorSet).neuralNet;
            save(fullfile(folderName, sensorID+"-linear-"+predictorSet+"-"+modelSaveDate+".mat"), "linearModel");
            save(fullfile(folderName, sensorID+"-net-"+predictorSet+"-"+modelSaveDate+".mat"),    "neuralNetModel");
        end
    end
end

function [alignedRefData, alignedDaqData] = alignDataWithCrossCorr(refData, daqData)

    if ~istimetable(refData) || ~istimetable(daqData)
        error('Both refData and daqData must be timetables.');
    end
    
    if ~ismember('C', refData.Properties.VariableNames)
        error('refData must contain a variable named ''C''.');
    end
    if ~ismember('CA', daqData.Properties.VariableNames)
        error('daqData must contain a variable named ''CA''.');
    end
    
    refTimeDiffs = seconds(diff(refData.T));
    daqTimeDiffs = seconds(diff(daqData.T));
    
    refMedianDT = median(refTimeDiffs);
    daqMedianDT = median(daqTimeDiffs);
    
    commonDT = min([refMedianDT, daqMedianDT]);
    refData = retime(refData, 'regular', 'linear', 'TimeStep', seconds(commonDT));
    daqData = retime(daqData, 'regular', 'linear', 'TimeStep', seconds(commonDT));
    
    startT = max([min(daqData.T), min(refData.T)]);
    endT = min([max(daqData.T), max(refData.T)]);

    daqData = daqData(daqData.T > startT & daqData.T < endT, :);
    refData = refData(refData.T > startT & refData.T < endT, :);

    refVec = refData.C;
    daqVec = daqData.CA;
    
    validIdx = ~isnan(refVec) & ~isnan(daqVec);
    refVec_clean = refVec(validIdx);
    daqVec_clean = daqVec(validIdx);
    
    refVec_clean = refVec_clean(:);
    daqVec_clean = daqVec_clean(:);
    
    [c, lags] = xcorr(daqVec_clean, refVec_clean, 'coeff');
    
    [~, idxMax] = max(abs(c));
    bestLag = lags(idxMax);
    
    timeShiftSeconds = bestLag * commonDT;
    
    if bestLag > 0
        shiftedDaqTime = daqData.T + seconds(timeShiftSeconds);
    elseif bestLag < 0
        shiftedDaqTime = daqData.T + seconds(timeShiftSeconds);
    else
        shiftedDaqTime = daqData.T;
    end
    
    shiftedDaqData = daqData;
    shiftedDaqData.T = shiftedDaqTime;
    
    syncData = synchronize(refData, shiftedDaqData, 'intersection', 'linear');
    syncedRef = syncData(:, "C");
    syncedDaq = syncData;
    syncedDaq.C = [];

    
    validRows = ~isnan(syncedRef.C) & ~isnan(syncedDaq.CA);
    alignedRefData = syncedRef(validRows, :);
    alignedDaqData = syncedDaq(validRows, :);
    
    fprintf('Best Lag: %d samples\n', bestLag);
    fprintf('Time Shift Applied to DAQ Data: %.2f seconds\n', timeShiftSeconds);

    fig = figure();
    subplot(2, 1, 1);
    plot(daqData.T, daqData.CA, 'r.',"DisplayName", "DAQ");
    hold on;
    plot(refData.T, refData.C, 'b.',"DisplayName", "Reference");
    ylabel("CO_2 (ppm)")
    xlabel("Time")
    title("Pre-Adjustment")
    legend('Location','best')

    subplot(2, 1, 2);
    plot(alignedDaqData.T, alignedDaqData.CA, 'r.',"DisplayName", "DAQ");
    hold on;
    plot(alignedRefData.T, alignedRefData.C, 'b.',"DisplayName", "Reference");
    ylabel("CO_2 (ppm)")
    xlabel("Time")
    title("Post-Adjustment")
    legend('Location','best')

    filename = sprintf("figs_training/CORRECTION_TIMESTAMP_%s_TO_%s.fig", string(startT, "MM_dd_yy"), string(endT, "MM_dd_yy"));
    saveas(fig, filename);
    close all;

end

function evaluateSavedModels(evaluationData, predictors, targetVariable)
    if istimetable(evaluationData)
        timeVec = evaluationData.Time; 
        dataTable = timetable2table(evaluationData);
    else
        if ~ismember('Time', evaluationData.Properties.VariableNames)
            error('If passing a table, it must contain a "Time" variable for plotting.');
        end
        timeVec = evaluationData.Time;
        dataTable = evaluationData;
    end

    [modelFiles, modelPath] = uigetfile('*.mat', ...
        'Select one or more model files to evaluate', ...
        'MultiSelect','on');

    if isequal(modelFiles,0)
        disp('No model files selected. Evaluation canceled.');
        return;
    end

    if ischar(modelFiles) || isstring(modelFiles)
        modelFiles = {modelFiles}; 
    end

    fig = figure('Name','Model Evaluation','Position',[200 200 1600 800]);
    hold on; grid on;
    plot(timeVec, dataTable.(targetVariable), 'k.', 'DisplayName','True Target');
    
    legendEntries = {'True Target'};

    for f = 1:numel(modelFiles)
        fullPath = fullfile(modelPath, modelFiles{f});
        loadedStruct = load(fullPath);

        if isfield(loadedStruct, 'linearModel')
            modelObj = loadedStruct.linearModel;
            modelType = 'Linear';
        elseif isfield(loadedStruct, 'neuralNetModel')
            modelObj = loadedStruct.neuralNetModel;
            modelType = 'Neural Net';
        else
            disp(['Skipping file (not recognized as a linear or net model): ', modelFiles{f}]);
            continue;
        end

        try
            switch modelType
                case 'Linear'
                    predictedVals = predict(modelObj, dataTable{:, predictors});
                case 'Neural Net'
                    predictedVals = modelObj(dataTable{:, predictors}');
                    predictedVals = predictedVals(:);
            end
        catch ME
            disp(['Error while predicting with model in file: ' modelFiles{f}]);
            disp(['Message: ' ME.message]);
            continue;
        end

        trueVals = dataTable.(targetVariable);
        rmseVal = sqrt(mean((predictedVals - trueVals).^2, 'omitnan'));
        ssRes   = sum((trueVals - predictedVals).^2, 'omitnan');
        ssTot   = sum((trueVals - mean(trueVals,'omitnan')).^2, 'omitnan');
        r2Val   = 1 - (ssRes/ssTot);

        fprintf('\nModel File: %s\n', modelFiles{f});
        fprintf('Model Type: %s\n', modelType);
        fprintf('RMSE: %.3f\n', rmseVal);
        fprintf('R^2:  %.3f\n', r2Val);

        plot(timeVec, predictedVals, '.', 'DisplayName', sprintf('%s (%s)\nRMSE=%.2f,R^2=%.2f', ...
            modelFiles{f}, modelType, rmseVal, r2Val));

        legendEntries{end+1} = sprintf('%s', modelFiles{f});
    end

    legend(legendEntries, 'Interpreter', 'none', 'Location', 'best');
    xlabel('Time'); ylabel(targetVariable);
    title('Model Evaluations on Provided Dataset');
end
