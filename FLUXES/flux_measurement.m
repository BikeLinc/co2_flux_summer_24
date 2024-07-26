% This script is for applying flux calculations to datasets for static and
% dynamic chambers.

% July 2024
%

% Clear workspace
clc, clear, close all

% Adding more functions to path from utility folder.
addpath(genpath('../UTILS'));

% Adding analysis configuration variables, for consistency across
% calculations.
cfg = config();

%% Date Range
%  Define the dates you want to look at fluxes for.

% Dates (Y, M, D) that you want to look at fluxes for.
dateStart = datetime(2024, 7, 8);
dateEnd = datetime(2024, 7, 20);

% Times that you want to looka t fluxes for.
start = timeofday(datetime("00:00", 'InputFormat', 'HH:mm')) + dateStart;
stop = timeofday(datetime("11:59", 'InputFormat', 'HH:mm')) + dateEnd;

% Mark true what flux chamber datasets you want to look at during the start
% ann stop duration range.
static = true;
dynamic = true;

%% Smooth Values
% These values dictate the delta between datapoints and the smoothing
% window that is used across the dataset.

if static
    licorRetime = seconds(30);
    licorWindow = minutes(5);
end

if dynamic
    daqRetime = seconds(30);
    daqWindow = minutes(5);
end

%% Collect Data
% Looks through flux datasets for ones matching datetimes.

% import metadata 
meta = readtable("../DATA/FLUX/listing.xlsx");

licorSets = {};
daqSets = {};

% for all days in range
for day = dateStart:dateEnd

    % get folder based on day
    folder = sprintf('%02d',day.Month) + "." + sprintf('%02d',day.Day) + "." + day.Year;
    
    % get files in folder
    directory = dir("../DATA/FLUX/" + folder);

    % for all files in folder
    for i = 3:numel(directory)
        file = directory(i);

        % select meta data row for current file-folder combination
        row = strcmp(meta{:,1}, file.name) & strcmp(meta{:,2}, folder);
        fileMeta = meta(row, :);
        
        % import licor file
        if file.name == "licor.txt"
    
            licor = IMPORTLICORFILE(file.folder + "/" + file.name);
            licor.META_ID = ones(height(licor),1)*find(row);
            licor.C_CALIB = licor.C*0.9883+24.6914;
            
            time = licor.T < stop - minutes(3) & licor.T > start - minutes(3);
            licor = licor(time, :);

            licorSets = [licorSets; {licor}];
        
        % import DAQ-like file
        elseif file.name(end-3:end) == ".txt" || file.name(end-3:end) == ".TXT"

            daq = IMPORTDAQFILE(file.folder + "/" + file.name);
            daq.META_ID = ones(height(daq),1)*find(row);
            
            time = daq.T < stop & daq.T > start;
            daq = daq(time, :);

            if height(daq) ~= 0
                daqSets = [daqSets; {daq}];
            end

        end
    end

end

% remove variables if turned off
if ~static
    licorSets = [];
end

if ~dynamic
    daqSets = [];
end

if isempty(daq)
    dynamic = false;
end

if isempty(licor)
    static = false;
end

%% Plot Raw Data
figure();
subplot(2, 2, 1);
hold on;
if static
    for i = 1:height(licorSets)
        licor = licorSets{i};
        plot(licor.T, licor.C, '.', 'DisplayName', "LICOR Raw - File " + mean(licor.META_ID));
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend();
title("Raw Static Chamber Data")
subplot(2, 2, 2);
hold on;



if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
        plot(daq.T, daq.CA, '.', 'DisplayName', "DAQ CA Raw - File " + mean(daq.META_ID));
        plot(daq.T, daq.CB, '.', 'DisplayName', "DAQ CB Raw - File " + mean(daq.META_ID));
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend();
title("Raw Dynamic Chamber Data")


%% smooth data

if static
    licorRaw= licor;
    licor = retime(licor,"regular", 'mean', 'TimeStep', licorRetime);
    licor = smoothdata(licor, 'movmean', licorWindow);
    licor = unique(licor);
end

if dynamic
    daqRaw = daq;
    daq = retime(daq,"regular", 'mean', 'TimeStep', daqRetime);
    daq = smoothdata(daq, 'movmean', daqWindow);
    daq = unique(daq);
end


%% Apply Calibrations

if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
        row = mean(daq.META_ID);
        metaData = meta(row,:);
           
        if (metaData.PAD_A ~= "UPDATE" || metaData.PAD_B ~= "UPDATE")

            modelCA = load("../ELT CALIB/CALIB_W_LICOR/models/C-CL-linear-15-Jul-2024.mat");
            modelCB = load("../ELT CALIB/CALIB_W_LICOR/models/E-CL-linear-15-Jul-2024.mat");
            
            % Linear
            daq.CA_CALIB = predict(modelCA.model_1_lin, [daq.CA]);
            daq.CB_CALIB = predict(modelCB.model_2_lin, [daq.CB]);
            
            % Network
            %daq.CA_CALIB = modelCA.model_1_net(daq.CA')';
            %daq.CB_CALIB = modelCB.model_2_net(daq.CB')';
    
        else
    
            daq.CA_CALIB = daq.CA;
            daq.CB_CALIB = daq.CB;
        end
    end
end

%% Plot Smoothed and Calibrated Data

subplot(2, 2, 3);
hold on;
if static
    for i = 1:height(licorSets)
        licor = licorSets{i};
        plot(licor.T, licor.C_CALIB,'.', 'DisplayName', "LICOR - File " + mean(licor.META_ID));
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend();
title("Processed Static Chamber Data")
subplot(2, 2, 4);
hold on;
if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
        plot(daq.T, daq.CA_CALIB,'.', 'DisplayName', "DAQ CA - File " + mean(daq.META_ID));
        plot(daq.T, daq.CB_CALIB,'.', 'DisplayName', "DAQ CB - File " + mean(daq.META_ID));
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend();
title("Processed Dynamic Chamber Data")
sgtitle("Raw and Processed Data")

%% Calculate Static Flux

if static
    licor.DCDT = [nan;diff(licor.C_CALIB)];
    licor = rmmissing(licor);
    
    licor.DCDT_mg = cfg.ppm_to_mg(licor.DCDT);
    licor.C_mg = cfg.ppm_to_mg(licor.C_CALIB);
    licor.F_mg = (licor.DCDT_mg.*cfg.V)./cfg.As;
end


%% Calculate Dynamic Flux

if dynamic

    Q = cfg.lpm_to_cms(daq.Q/1000);% mLpm to cms
    AMBIENT = cfg.ppm_to_mg(daq.CB_CALIB);
    CHAMBER = cfg.ppm_to_mg(daq.CA_CALIB);
    CHAMBER_FLOOR = CHAMBER - AMBIENT(1);
    
    DCDT = [nan;diff(CHAMBER)];
    V = cfg.V;
    As = cfg.As;
    
    daq.F_mg = ((V.*DCDT) + Q.*(CHAMBER-AMBIENT))./As;
    
    DCDT_FLOOR = [nan;diff(CHAMBER_FLOOR)];
    daq.F_FLOOR_mg = ((V.*DCDT) + Q.*(CHAMBER_FLOOR-AMBIENT))./As;
end

%% Expected Range

range = [0.5e-6 10e-6];% umol
range = cfg.mol_to_ppm(range);
range = cfg.ppm_to_mg(range);

%% Plot

figure();
hold on;

if static
    plot(licor.T, licor.F_mg,'r^', 'DisplayName', "LICOR");
    yline(mean(licor.F_mg),'r--', 'DisplayName', "LICOR \mu = " + mean(licor.F_mg) + " mg/m^2/s")
end

if dynamic
    plot(daq.T, daq.F_mg,'gs', 'DisplayName', "DAQ");
    plot(daq.T, daq.F_FLOOR_mg,'b*', 'DisplayName', "DAQ Floor");
    yline(mean(rmmissing(daq.F_mg)),'g--', 'DisplayName', "DAQ \mu = " + mean(rmmissing(daq.F_mg)) + " mg/m^2/s")
    yline(mean(rmmissing(daq.F_FLOOR_mg)),'b--', 'DisplayName', "DAQ Floor \mu = " + mean(rmmissing(daq.F_FLOOR_mg)) + " mg/m^2/s")
end

xlabel("Time");
ylabel("Flux Rate [mg/m^2/s]");
title("Measured CO_2 Flux Rates  (" + start.Month + "/" + start.Day + "/" + start.Year + " to " + stop.Month + "/" + stop.Day + "/" + stop.Year + ")");
grid on;
legend('location','eastoutside');
yregion(range, 'FaceColor',"magenta",'FaceAlpha', 0.1, 'DisplayName', 'Expected Range');

%% Generate Report

clc
fprintf("---- Flux Measurement Report ----\n")

fprintf("\nTime Range\n\n")
fprintf("\tStart:\t%s\n",start)
fprintf("\tStop:\t%s\n",stop)

if static
    fprintf("\nStatic Chamber LICOR\n\n")
    fprintf("\tNo. Points:\t%d\n", numel(licor))
    fprintf("\tRetimed:\t%s\n",licorRetime);
    fprintf("\tWindow:\t\t%s\n",licorWindow);
    
    fprintf("\nStatic Chamber DAQ\n\n")
    fprintf("\tNo. Points:\t%d\n", numel(daq))
    fprintf("\tRetimed:\t%s\n",daqRetime);
    fprintf("\tWindow:\t\t%s\n",daqWindow);
    
    fprintf("\nStatic Calibrations\n\n")
    fprintf("\tLICOR:\t0.9883x+24.6914 ppm\n")
    
    
    fprintf("\nStatic Fluxes\n\n")
    fprintf("\tStatic:\t\t\t%03f mg/m^2/s\n", mean(licor.F_mg))
end


if  dynamic
    fprintf("\tCB:\t\tLinear R^2 %0.4f\n", modelCB.model_2_lin.Rsquared.Ordinary)
    fprintf("\tCA:\t\tLinear R^2 %0.4f\n", modelCA.model_1_lin.Rsquared.Ordinary)
    fprintf("\tDynamic:\t\t%03f mg/m^2/s\n", mean(rmmissing(daq.F_mg)))
    fprintf("\tDynamic Floor:\t%03f mg/m^2/s\n", mean(rmmissing(daq.F_FLOOR_mg)))
end