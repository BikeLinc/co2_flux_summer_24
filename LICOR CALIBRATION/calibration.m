%% Generate A Calibration for LICOR Dataset

licor = IMPORTLICORFILE("calibration_dataset_licor.data");

%%
t_400_start = datetime("today") + hours(11) + minutes(4);
t_400_end = datetime("today") + hours(11) + minutes(5);
t_3003_start = datetime("today") + hours(11) + minutes(10);
t_3003_end = datetime("today") + hours(11) + minutes(11);
writetimetable(licor)


%%

data = readtable("licor.csv");
lm = fitlm(data.C, data.CT);
plot(lm)
