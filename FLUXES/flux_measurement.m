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
static = true;
dynamic = true;

%% Data Settings
% These values dictate the delta between datapoints and the smoothing
% window that is used across the dataset.

if static
    licorRetime = seconds(5);
    licorWindow = minutes(1);
end

if dynamic
    daqRetime = seconds(500);
    daqWindow = minutes(1);
end

%% Collect Data
% Looks through flux datasets for files matching datetime range. Looks
% through listing.xlsx to collect metadata for deployment information,
% sensor location, and calibrations required.

% import metadata
meta = readtable("../DATA/FLUX/listing.xlsx");

% variables for collecting data and metadata
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


            licor.data = IMPORTLICORFILE(file.folder + "/" + file.name);     % call import function
            metaID = ones(height(licor),1)*find(row);                        % get metadata id
            licor.data.C_CALIB = licor.data.C*0.9883+24.6914;                % apply licor calibration
            licor.meta = meta(metaID, :);                                    % link metadata


            % confine datetime to desired date range
            time = licor.data.T < stop - minutes(3) & licor.data.T > start - minutes(3);
            licor.data = licor.data(time, :);
            licor.file = file.folder + "/" + file.name;

            licorSets = [licorSets; {licor}];

            % import DAQ-like file
        elseif file.name(end-3:end) == ".txt" || file.name(end-3:end) == ".TXT"

            daq.data = IMPORTDAQFILE(file.folder + "/" + file.name);
            metaID = ones(height(daq),1)*find(row);

            time = daq.data.T < stop & daq.data.T > start;
            daq.data = daq.data(time, :);

            daq.meta = meta(metaID, :);                                    % link metadata

            daq.file = file.folder + "/" + file.name;

            disp(sprintf(file.name + "\t" + height(daq.data)));

            if height(daq.data) ~= 0
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

if isempty(daqSets)
    dynamic = false;
end

if isempty(licorSets)
    static = false;
end

% count number of datasets
licorCount = height(licorSets);
daqCount = height(daqSets);





%% Plot Raw Data

figure();
subplot(3, 2, 1);
hold on;

h0 = [];
for i = 1:height(licorSets)
    licor = licorSets{i};
    p = plot(licor.data.T, licor.data.C, 'b.');
    h0 = [h0, p];
end

ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend(h0, "Raw LICOR 7810");
title("Raw Static Chamber Data")
subplot(3, 2, 2);
hold on;


h0 = [];
h1 = [];
if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
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

if static
    for i = 1:height(licorSets)
        licor = licorSets{i};
        licorRaw= licor.data;
        licor.data = retime(licor.data,"regular", 'mean', 'TimeStep', licorRetime);
        licor.data = smoothdata(licor.data, 'movmean', licorWindow);
        licor.data = unique(licor.data);
        licorSets{i} = licor;
    end


end

if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
        daqRaw = daq.data;
        daq.data = retime(daq.data,"regular", 'mean', 'TimeStep', daqRetime);
        daq.data = smoothdata(daq.data, 'movmean', daqWindow);
        daq.data = unique(daq.data);
        daqSets{i} = daq;

    end

end

%% Calculate Ambient Sensor Automatically
% Beacuse in the field we can occationally not document which sensor is
% attached to what, we use a simple standard deviation after the smoothing
% process to determine what sensor was the ambient and what sensor was the
% chamber. We then flip the datasets accordingly. Ambient sensor should
% have a lower STD that the chamber because the chamber must reach steady
% state after deployment.

if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};

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

        daqSets{i} = daq;

    end
end

%% Plot Smoothed & Auto-Corr Data

subplot(3, 2, 3);
hold on;

h0 = [];
for i = 1:height(licorSets)
    licor = licorSets{i};
    p = plot(licor.data.T, licor.data.C, 'b.');
    h0 = [h0, p];
end

ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend(h0, "Smoothed LICOR 7810");
title("Smoothed Static Chamber Data")
subplot(3, 2, 4);
hold on;


h0 = [];
h1 = [];
if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
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

if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};

        disp(daq.meta.PAD_A);

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

        daqSets{i} = daq;
    end
end

%% Plot Calibrated Data

subplot(3, 2, 5);
hold on;

h0 = [];
for i = 1:height(licorSets)
    licor = licorSets{i};
    p = plot(licor.data.T, licor.data.C_CALIB, 'b.');
    h0 = [h0, p];
end

ylabel("CO_2 [ppm]")
xlabel("Time")
grid on;
legend(h0, "Calibrated LICOR 7810");
title("Calibrated Static Chamber Data")
subplot(3, 2, 6);
hold on;


h0 = [];
h1 = [];
if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
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

if static
    for i = 1:height(licorSets)
        licor = licorSets{i};
        dc = [0; diff(licor.data.C_CALIB)];
        dt = [0; seconds(diff(licor.data.T))];

        licor.data.DCDT = dc./dt;
        licor.data = rmmissing(licor.data);

        licor.data.DCDT_mg = cfg.ppm_to_mg(licor.data.DCDT);
        licor.data.C_mg = cfg.ppm_to_mg(licor.data.C_CALIB);
        licor.data.F_mg = (licor.data.DCDT_mg.*cfg.V)./cfg.As;
        licorSets{i} = licor;
    end


end


%% Calculate Dynamic Flux

if dynamic

    for i = 1:height(daqSets)
        daq = daqSets{i};

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

        daqSets{i} = daq;
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
if static
    for i = 1:height(licorSets)
        licor = licorSets{i};
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
if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};
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


licorBurnLump = [];
licorUnburnLump = [];
daqUnburnLump = [];
daqBurnLump = [];

if static
    for i = 1:height(licorSets)
        licor = licorSets{i};
        if (licor.meta.DEPLOY_SITE == "BURN")
            licorBurnLump = [licorBurnLump; licor.data];
        else
            licorUnburnLump = [licorUnburnLump; licor.data];
        end
    end
end

if dynamic
    for i = 1:height(daqSets)
        daq = daqSets{i};

        if (daq.meta.DEPLOY_SITE == "BURN")
            daqBurnLump = [daqBurnLump; daq.data];
        else
            daqUnburnLump = [daqUnburnLump; daq.data];
        end

    end
end

%% Measured Fluxes Over a Day Long Period
daqBurnLumpDay=daqBurnLump.F_mg*86400;
daqUnburnLumpDay=daqUnburnLump.F_mg*86400;
licorBurnLumpDay=licorBurnLump.F_mg*86400;
licorUnburnLumpDay=licorUnburnLump.F_mg*86400;
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
disp("Burned DAQ Mean " + mean(rmmissing(daqBurnLumpDay)))
disp("Unburned DAQ Mean " + mean(rmmissing(daqUnburnLumpDay)))
disp("Burned Licor Mean " + mean(rmmissing(licorBurnLumpDay)))
disp("Unburned Licor Mean " + mean(rmmissing(licorUnburnLumpDay)))

disp("Burned DAQ Median " + median(rmmissing(daqBurnLumpDay)))
disp("Unburned DAQ Median " + median(rmmissing(daqUnburnLumpDay)))
disp("Burned Licor Median " + median(rmmissing(licorBurnLumpDay)))
disp("Unburned Licor Median " + median(rmmissing(licorUnburnLumpDay)))

%% Generate Report


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
    fprintf("\tStatic:\t\t\t%03f mg/m^2/s\n", mean(licor.data.F_mg))
end


if  dynamic
    %fprintf("\tCB:\t\tLinear R^2 %0.4f\n", modelCB.model_2_lin.Rsquared.Ordinary)
    %fprintf("\tCA:\t\tLinear R^2 %0.4f\n", modelCA.model_1_lin.Rsquared.Ordinary)
    %fprintf("\tDynamic:\t\t%03f mg/m^2/s\n", mean(rmmissing(daq.F_mg)))
    %fprintf("\tDynamic Floor:\t%03f mg/m^2/s\n", mean(rmmissing(daq.F_FLOOR_mg)))
end

