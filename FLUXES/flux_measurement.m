% This script is for applying flux calculations to datasets for static and
% dynamic chambers. Sections generally use processDynamic and processStatic
% logical flags to run or not run specific sections.
%
% July 2024
%

% Clear workspace
clc, clear, close all

fprintf("Running:\tflux_measurement.m\n")

% Adding more functions to path from utility folder.
addpath(genpath('../UTILS'));

% Adding analysis configuration variables, for consistency across
% calculations.
cfg = config();

%% Date Range Settings
%  Define the dates you want to look at fluxes for.

% Dates (Y, M, D) that you want to look at fluxes for.
dateStart = datetime(2024, 7, 8);
dateEnd = datetime(2024, 7, 31);

% Times that you want to looka t fluxes for.
start = timeofday(datetime("0:59", 'InputFormat', 'HH:mm')) + dateStart;
stop = timeofday(datetime("11:59", 'InputFormat', 'HH:mm')) + dateEnd;

% Mark true what flux chamber datasets you want to look at during the start
% ann stop duration range.
processStatic = true;
processDynamic = true;

%% Data Settings
% These values dictate the delta between datapoints and the smoothing
% window that is used across the dataset.

if processStatic
    licorRetime = seconds(30);
    licorWindow = minutes(5);
end

if processDynamic
    daqRetime = seconds(30);
    daqWindow = minutes(60);
end

%% Collect Data
% Looks through flux datasets for files matching datetime range. Looks
% through listing.xlsx to collect metadata for deployment information,
% sensor location, and calibrations required.

% import metadata from table
metadata = readtable("../DATA/FLUX/listing.xlsx");

% variables for collecting data and metadata
licorDatasets = {};
daqDatasets = {};

fprintf("Importing Data:\n")

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
        row = strcmp(metadata{:,1}, file.name) & strcmp(metadata{:,2}, folder);
        fileMeta = metadata(row, :);

        % import licor file
        if file.name == "licor.txt"


            licor.data = IMPORTLICORFILE(file.folder + "/" + file.name);     % call import function
            metaID = ones(height(licor),1)*find(row);                        % get metadata id
            licor.data.C_CALIB = licor.data.C*0.9883+24.6914;                % apply licor calibration
            licor.meta = metadata(metaID, :);                                    % link metadata


            % confine datetime to desired date range
            time = licor.data.T < stop - minutes(3) & licor.data.T > start - minutes(3);
            licor.data = licor.data(time, :);
            licor.data = rmmissing(licor.data);
            licor.file = file.folder + "/" + file.name;

            fprintf("\tFile Loaded: (LICOR)\t%s/%s\n", folder, file.name);

            licorDatasets = [licorDatasets; {licor}];

            % import DAQ-like file
        elseif file.name(end-3:end) == ".txt" || file.name(end-3:end) == ".TXT"

            daq.data = IMPORTDAQFILE(file.folder + "/" + file.name);
            metaID = ones(height(daq),1)*find(row);

            time = daq.data.T < stop & daq.data.T > start;
            daq.data = daq.data(time, :);
            daq.data = rmmissing(daq.data);

            daq.meta = metadata(metaID, :);                                    % link metadata

            daq.file = file.folder + "/" + file.name;

            fprintf("\tFile Loaded: (DAQ)\t%s/%s\n", folder, file.name);

            if height(daq.data) ~= 0
                daqDatasets = [daqDatasets; {daq}];
            end

        end
    end

end

% remove variables if turned off
if ~processStatic
    licorDatasets = [];
end

if ~processDynamic
    daqDatasets = [];
end

if isempty(daqDatasets)
    processDynamic = false;
end

if isempty(licorDatasets)
    processStatic = false;
end

% count number of datasets
licorCount = height(licorDatasets);
daqCount = height(daqDatasets);





%% Plot Raw Data

figure();
subplot(3, 2, 1);
hold on;

h0 = [];
for i = 1:height(licorDatasets)
    licor = licorDatasets{i};
    p = plot(licor.data.T, licor.data.C, 'b.');
    h0 = [h0, p];
end

ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend(h0(1), "Raw LICOR 7810");
title("Raw Static Chamber Data")
subplot(3, 2, 2);
hold on;


h0 = [];
h1 = [];
if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};
        p0 = plot(daq.data.T, daq.data.CA, 'g.');
        p1 = plot(daq.data.T, daq.data.CB, 'm.');
        h0 = [h0, p0];
        h1 = [h1, p1];
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend([h0(1), h1(1)], ["Raw Pad A Sensor", "Raw Pad B Sensor"]);
title("Raw Dynamic Chamber Data")

%% Smooth Data

if processStatic
    for i = 1:height(licorDatasets)
        licor = licorDatasets{i};
        licorRaw= licor.data;
        licor.data = retime(licor.data,"regular", 'mean', 'TimeStep', licorRetime);
        licor.data = smoothdata(licor.data, 'movmean', licorWindow);
        licor.data = unique(licor.data);
        licorDatasets{i} = licor;
    end


end

if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};
        daqRaw = daq.data;
        daq.data = retime(daq.data,"regular", 'mean', 'TimeStep', daqRetime);
        daq.data = smoothdata(daq.data, 'movmean', daqWindow);
        daq.data = unique(daq.data);
        daqDatasets{i} = daq;

    end

end

%% Calculate Ambient Sensor Automatically
% Beacuse in the field we can occationally not document which sensor is
% attached to what, we use a simple standard deviation after the smoothing
% process to determine what sensor was the ambient and what sensor was the
% chamber. We then flip the datasets accordingly. Ambient sensor should
% have a lower STD that the chamber because the chamber must reach steady
% state after deployment.

if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};

        caSTD = std(daq.data.CA);
        cbSTD = std(daq.data.CB);

        if caSTD > cbSTD

            daqCB = daq.data.CB;
            daqTB = daq.data.TB;
            daqHB = daq.data.HB;

            daqCA = daq.data.CA;
            daqTA = daq.data.TA;
            daqHA = daq.data.HA;

            daq.data.CB = daqCA;
            daq.data.TB = daqTA;
            daq.data.HB = daqHA;

            daq.data.CA = daqCB;
            daq.data.TA = daqTB;
            daq.data.HA = daqHB;

        end

        daqDatasets{i} = daq;

    end
end

%% Plot Smoothed & Auto-Corr Data

subplot(3, 2, 3);
hold on;

h0 = [];
for i = 1:height(licorDatasets)
    licor = licorDatasets{i};
    p = plot(licor.data.T, licor.data.C, 'b.');
    h0 = [h0, p];
end

ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend(h0(1), "Smoothed LICOR 7810");
title("Smoothed Static Chamber Data")
subplot(3, 2, 4);
hold on;


h0 = [];
h1 = [];
if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};
        p0 = plot(daq.data.T, daq.data.CA, 'g.');
        p1 = plot(daq.data.T, daq.data.CB, 'm.');
        h0 = [h0, p0];
        h1 = [h1, p1];
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend([h0(1), h1(1)], ["Smoothed Auto-Corr Pad A Sensor", "Smoothed Auto-Corr Pad B Sensor"]);
title("Smoothed Dynamic Chamber Data")


%% Apply Calibrations

if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};

        if (daq.meta.PAD_A_ELT_SENSOR == "E" || daq.meta.PAD_B_ELT_SENSOR == "C")

            modelCA = load("../ELT CALIB/CALIB_W_LICOR/models/E-CL-linear-15-Jul-2024.mat");
            modelCB = load("../ELT CALIB/CALIB_W_LICOR/models/C-CL-linear-15-Jul-2024.mat");

            % Linear
            daq.data.CA_CALIB = predict(modelCA.model_2_lin, [daq.data.CA]);
            daq.data.CB_CALIB = predict(modelCB.model_1_lin, [daq.data.CB]);

            % Network
            %daq.CA_CALIB = modelCA.model_2_net(daq.CA')';
            %daq.CB_CALIB = modelCB.model_1_net(daq.CB')';

        else

            daq.data.CA_CALIB = daq.data.CA;
            daq.data.CB_CALIB = daq.data.CB;
        end

        daqDatasets{i} = daq;
    end
end

%% Plot Calibrated Data

subplot(3, 2, 5);
hold on;

h0 = [];
for i = 1:height(licorDatasets)
    licor = licorDatasets{i};
    p = plot(licor.data.T, licor.data.C_CALIB, 'b.');
    h0 = [h0, p];
end

ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend(h0(1), "Calibrated LICOR 7810");
title("Calibrated Static Chamber Data")
subplot(3, 2, 6);
hold on;


h0 = [];
h1 = [];
if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};
        p0 = plot(daq.data.T, daq.data.CA_CALIB, 'g.');
        p1 = plot(daq.data.T, daq.data.CB_CALIB, 'm.');
        h0 = [h0, p0];
        h1 = [h1, p1];
    end
end
ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend([h0(1), h1(1)], ["Calibrated Pad A Sensor", "Calibrated Pad B Sensor"]);
title("Calibrated Dynamic Chamber Data")



%% Calculate Static Flux

if processStatic
    for i = 1:height(licorDatasets)
        licor = licorDatasets{i};
        dc = [0; diff(licor.data.C_CALIB)];
        dt = [0; seconds(diff(licor.data.T))];

        licor.data.DCDT = dc./dt;
        licor.data = rmmissing(licor.data);

        licor.data.DCDT_mg = cfg.ppm_to_mg(licor.data.DCDT);
        licor.data.C_mg = cfg.ppm_to_mg(licor.data.C_CALIB);
        licor.data.F_mg = (licor.data.DCDT_mg.*cfg.V)./cfg.As;
        licorDatasets{i} = licor;
    end


end


%% Calculate Dynamic Flux

if processDynamic

    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};

        Q = cfg.lpm_to_cms(mean(daq.data.Q/1000));

        if Q == 0
            Q = cfg.lpm_to_cms(2.0204125);
        end

        AMBIENT = cfg.ppm_to_mg(daq.data.CA_CALIB);
        CHAMBER = cfg.ppm_to_mg(daq.data.CB_CALIB);
        CHAMBER_FLOOR = CHAMBER - AMBIENT(1);

        V = cfg.V;
        As = pi*(0.1603375^2);

        daq.data.F_mg = (Q.*(CHAMBER-AMBIENT))./As;

        daqDatasets{i} = daq;
    end
end

%% Expected Range

range = [0.5e-6 10e-6];% umol
range = cfg.mol_to_ppm(range);
range = cfg.ppm_to_mg(range);

%% Plot

figure();
hold on;

h0 = [];
h1 = [];
if processStatic
    for i = 1:height(licorDatasets)
        licor = licorDatasets{i};
        if (licor.meta.DEPLOY_SITE == "BURN")
            p0 = plot(licor.data.T, licor.data.F_mg,'r^');
            h0 = [h0; p0];
        else
            p1 = plot(licor.data.T, licor.data.F_mg,'m^');
            h1 = [h1; p1];
        end
    end
end

h2 = [];
h3 = [];
if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};
        if (daq.meta.DEPLOY_SITE == "BURN")
            p0 = plot(daq.data.T, daq.data.F_mg,'bs');
            h2 = [h2; p0];
        else
            p1 = plot(daq.data.T, daq.data.F_mg,'gs');
            h3 = [h3; p1];
        end
    end
end

xlabel("Time");
ylabel("Flux Rate [mg/m^2/s]");
title("Measured CO_2 Flux Rates  (" + start.Month + "/" + start.Day + "/" + start.Year + " to " + stop.Month + "/" + stop.Day + "/" + stop.Year + ")");
grid on;
legend([h0(1), h1(1), h2(1), h3(1)], ["Static Burn","Static Un-Burn", "Dynamic Burn", "Dynamic Un-Burn"]);
yregion(range, 'FaceColor',"magenta",'FaceAlpha', 0.1, 'DisplayName', 'Expected Range');


%% Lump Data Together


licorBurn = [];
licorUnburn = [];
daqUnburn = [];
daqBurn = [];

% collect all licor data
if processStatic
    for i = 1:height(licorDatasets)
        licor = licorDatasets{i};
        if (licor.meta.DEPLOY_SITE == "BURN")
            licorBurn = [licorBurn; licor.data];
        else
            licorUnburn = [licorUnburn; licor.data];
        end
    end
end
licorUnburn = rmmissing(licorUnburn);
licorBurn = rmmissing(licorBurn);

if processDynamic
    for i = 1:height(daqDatasets)
        daq = daqDatasets{i};

        if (daq.meta.DEPLOY_SITE == "BURN")
            daqBurn = [daqBurn; daq.data];
        else
            daqUnburn = [daqUnburn; daq.data];
        end

    end
end
daqUnburn = rmmissing(daqUnburn);
daqBurn = rmmissing(daqBurn);


%% Measured Fluxes Over a Day Long Period
daqBurnLumpDay=daqBurn.F_mg*86400;
daqUnburnLumpDay=daqUnburn.F_mg*86400;
licorBurnLumpDay=licorBurn.F_mg*86400;
licorUnburnLumpDay=licorUnburn.F_mg*86400;
%% Generate Box-And-Whisker Plots

% create figure, hold on, and add grid
figure();
hold on;
grid on;

% lumped data to plot
x1 = daqBurnLumpDay;
x2 = daqUnburnLumpDay;
x3 = licorBurnLumpDay;
x4 = licorUnburnLumpDay;
x = [x1; x2; x3; x4];

% lumped categories to plot
g1 = repmat({'DAQ Burned'}, height(x1),1);
g2 = repmat({'DAQ Un-Burned'}, height(x2),1);
g3 = repmat({'LICOR Burned'}, height(x3),1);
g4 = repmat({'LICOR Un-Burned'}, height(x4),1);
g = [g1; g2; g3; g4];

% plot box and whisker
boxplot(x, g);
xlabel("Location Analyzed and Instrument Used");
ylabel("Flux Rate [mg/m^2/day]");
title("Summarized Data Distribution of Measured CO_2 Flux Rates");
yregion(range*86400, 'FaceColor',"magenta",'FaceAlpha', 0.1, 'DisplayName', 'Expected Range');



%% Generate Written Summary


fprintf("Reporting Data:\n")

fprintf("\tSITE\tSOURCE\tMEAN\tMEDIAN\t (mg/m^2/s)\n")
fprintf("\t%s\t%s\t%0.4f\t%0.4f\n","Burn", "DAQ", mean(daqBurnLumpDay), median(daqBurnLumpDay))


% disp("Burned DAQ Mean " + mean(rmmissing(daqBurnLumpDay)))
% disp("Unburned DAQ Mean " + mean(rmmissing(daqUnburnLumpDay)))
% disp("Burned Licor Mean " + mean(rmmissing(licorBurnLumpDay)))
% disp("Unburned Licor Mean " + mean(rmmissing(licorUnburnLumpDay)))
% disp("Burned DAQ Median " + median(rmmissing(daqBurnLumpDay)))
% disp("Unburned DAQ Median " + median(rmmissing(daqUnburnLumpDay)))
% disp("Burned Licor Median " + median(rmmissing(licorBurnLumpDay)))
% disp("Unburned Licor Median " + median(rmmissing(licorUnburnLumpDay)))