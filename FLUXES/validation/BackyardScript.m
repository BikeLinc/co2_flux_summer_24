%%
clc, clear, close all
%%

data = readtable("data/BackyardData/BackyardTestData.csv");
data = table2timetable(data);

%%

config = analysis_config();

%%

load("calib.mat");

data = retime(data, 'regular', 'mean', 'TimeStep', minutes(5));

%%
data.CA = predict(lin_rega, data.CA);
data.CB = predict(lin_regb, data.CB);


%%

figure();
hold on
%plot(data.T, data.CA, 'DisplayName', "CA")
%plot(data.T, data.CB, 'DisplayName', "CB")

Q = config.lpm_to_cms(1.227*(data.W/1000)+0.0143);
CA = config.ppm_to_mol(data.CA);
CB = config.ppm_to_mol(data.CB);

f = ((Q.*(CB-CA))./config.As)*1e6;

plot(data.T, f, 'DisplayName', "Flux")
legend();