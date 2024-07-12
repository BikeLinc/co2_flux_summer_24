clc, clear, close all

addpath(genpath('../UTILS'));

%% Import Data

licor = IMPORTLICORFILE("../DATA/FLUX/07.10.2024/licor_7_10_24.txt");
licor = rmmissing(licor);



%% smooth data
licor_sm = retime(licor,"regular", 'mean', 'TimeStep', seconds(30));
licor_sm = smoothdata(licor_sm, 'movmean', minutes(2));
licor_sm = unique(licor_sm);


%% Calculate Differences

licor_sm.DCDT = [nan;diff(licor_sm.C)];
licor_sm = rmmissing(licor_sm);

%% Plot
figure();
plot(licor_sm.T, licor_sm.DCDT, 'DisplayName', 'dC/dt');
hold on;
plot(licor_sm.T, licor_sm.C, 'DisplayName', 'C');