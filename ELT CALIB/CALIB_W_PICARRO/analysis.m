%% analysis.m
%
% A tool for analysing PICARRO and DAQ colocations.
%
% May 2024
%

clc, clear, close all
%% Import DAQ Dataset

daq = importdaqfile("daq_gh_colo2.txt");


%% Import PICARRO Datasets from subfolder


if isfile("picarro.txt")
picarro = readtimetable("picarro.txt");
else
    picarro = importpicarrosubfolders("picarro/2024");
end

%% Clean Datasets
daq = rmmissing(daq);
picarro = rmmissing(picarro);


%% Time Overlap and Shift

% Shift Picarro Dataset By 6 Hours
% Corresponds to GMT (-7) and Daylight Savings (+1) Time Shifts
picarro.T = picarro.T - hours(6);


%% Find Overlapping Times
daq_max_t = max(daq.T);
daq_min_t = min(daq.T);
picarro_max_t = max(picarro.T);
picarro_min_t = min(picarro.T);

% confine picarro dataset to daq range
picarro_idx = picarro.T < daq_max_t & picarro.T > daq_min_t;
picarro = picarro(picarro_idx, :);

%% Plot RAW Data
fig1  = figure();
hold on
plot(picarro.T, picarro.CO2_sync)
plot(daq.T, daq.CA)
plot(daq.T, daq.CB)
legend(["Picarro CO2 (Sync)", "DAQ ELT A", "DAQ ELT B"])
xlabel("Datetime")
ylabel("CO_2 [ppm]")
title("Ambient Colocation Raw Data")
grid on

%% Apply Windowing Filter and Syncronize Dataset
windowSize = 1000;
num_trans = (1/windowSize)*ones(1,windowSize);
den_trans = 1;
picarro.CO2_sync = filter(num_trans, den_trans, picarro.CO2_sync);
daq.CA = filter(num_trans, den_trans, daq.CA);
daq.CB = filter(num_trans, den_trans, daq.CB);
picarro(1:windowSize*2, :) = [];
daq(1:windowSize*2, :) = [];

data = synchronize(daq, picarro);
data = rmmissing(data);

%% Plot Smooth Data
fig2  = figure();
hold on
plot(data.T, data.CO2_sync)
plot(data.T, data.CA)
plot(data.T, data.CB)
legend(["Picarro CO2 (Sync)", "DAQ ELT A", "DAQ ELT B"])
xlabel("Datetime")
ylabel("CO_2 [ppm]")
title("Ambient Colocation Smoothed Data")
grid on


%% Split Data for Regressions
cv = cvpartition(size(data,1), 'HoldOut', 0.3);
data = [data.CA, data.TA, data.HA, data.CB, data.TB, data.HB, data.CO2_sync];
testd = data(cv.test, :);
traind = data(~cv.test, :);
%%
train_ca = traind(:, 1:3);
train_cb = traind(:, 4:6);
train_ct = traind(:, 7);
test_ca = testd(:, 1:3);
test_cb = testd(:, 4:6);
test_ct = testd(:, 7);

% %% Linear Regression
 lin_rega = fitlm(train_ca, train_ct);
 lin_regb = fitlm(train_cb, train_ct);
% 
% step = 100;
% results = [];
% best_lin = {};
% best_ann = [];
% best_r2_ann = [];
% best_r2_lin = {};
% best_lin_rmse = inf;
% best_ann_rmse = inf;
% best_ann_r2_avg = 0;
% best_lin_r2_avg = 0;
% best_lin_boundry = min(train_ct);
% best_ann_boundry = min(train_ct);
% best_r2_ann_boundry = min(train_ct);
% best_r2_lin_boundry = min(train_ct);
% 
% 
% for boundry = min(train_ct)+step:step:max(train_ct)
% 
%     disp(100*(boundry-min(train_ct))/max(train_ct) + "% Complete" )
% 
% 
%     if (length(train_ca(train_ct > boundry)) > 1)
% 
%         [lin_rega_low, lin_rega_low_gof, ~] = fit(train_ca(train_ct < boundry), train_ct(train_ct < boundry), 'poly2');
%         [lin_rega_high, lin_rega_high_gof, ~] = fit(train_ca(train_ct >= boundry), train_ct(train_ct >= boundry), 'poly2');
%         [lin_regb_low, lin_regb_low_gof, ~] = fit(train_cb(train_ct < boundry), train_ct(train_ct < boundry), 'poly2');
%         [lin_regb_high, lin_regb_high_gof, ~] = fit(train_cb(train_ct >= boundry), train_ct(train_ct >= boundry), 'poly2');
% 
%         lin_rega_high_pred = lin_rega_high(test_ca(test_ct > boundry));
%         lin_regb_high_pred = lin_rega_high(test_ca(test_ct > boundry));
%         lin_rega_low_pred = lin_rega_low(test_ca(test_ct < boundry));
%         lin_regb_low_pred = lin_rega_low(test_ca(test_ct < boundry));
% 
%         %disp("A Low:  " + lin_rega_low.RMSE)
%         %disp("A High: " + lin_rega_high.RMSE)
%         %disp("B Low:  " + lin_regb_low.RMSE)
%         %disp("B High: " + lin_regb_high.RMSE)
%         %disp("-")
% 
% 
%         ann_rega_low = feedforwardnet([16, 16]);
%         ann_regb_low = feedforwardnet([16, 16]);
%         ann_rega_high = feedforwardnet([16, 16]);
%         ann_regb_high = feedforwardnet([16, 16]);
%         ann_rega_low.trainParam.showWindow = false;
%         ann_regb_low.trainParam.showWindow = false;
%         ann_rega_high.trainParam.showWindow = false;
%         ann_regb_high.trainParam.showWindow = false;
%         ann_rega_low = train(ann_rega_low, train_ca(train_ct < boundry)', train_ct(train_ct < boundry)');
%         ann_regb_low = train(ann_regb_low, train_cb(train_ct < boundry)', train_ct(train_ct < boundry)');
%         ann_rega_high = train(ann_rega_high, train_ca(train_ct > boundry)', train_ct(train_ct > boundry)');
%         ann_regb_high = train(ann_regb_high, train_cb(train_ct > boundry)', train_ct(train_ct > boundry)');
% 
%         ann_rega_high_pred = ann_rega_high(test_ca(test_ct > boundry)')';
%         ann_regb_high_pred = ann_regb_high(test_cb(test_ct > boundry)')';
%         ann_rega_low_pred = ann_rega_low(test_ca(test_ct < boundry)')';
%         ann_regb_low_pred = ann_regb_low(test_cb(test_ct < boundry)')';
% 
%         ann_rega_low_rmse = sqrt(mean((ann_rega_low_pred - test_ct(test_ct < boundry)).^2));
%         ann_regb_low_rmse = sqrt(mean((ann_regb_low_pred - test_ct(test_ct < boundry)).^2));
%         ann_rega_high_rmse = sqrt(mean((ann_rega_high_pred - test_ct(test_ct > boundry)).^2));
%         ann_regb_high_rmse = sqrt(mean((ann_regb_high_pred - test_ct(test_ct > boundry)).^2));
% 
%         ann_rega_r2_low = 1 - ((sum((ann_rega_low_pred - test_ct(test_ct < boundry)).^2))/(sum(((test_ct(test_ct < boundry) - mean(test_ct(test_ct < boundry))).^2))));
%         ann_regb_r2_low = 1 - ((sum((ann_regb_low_pred - test_ct(test_ct < boundry)).^2))/(sum(((test_ct(test_ct < boundry) - mean(test_ct(test_ct < boundry))).^2))));
%         lin_rega_r2_low = 1 - ((sum((lin_rega_low_pred - test_ct(test_ct < boundry)).^2))/(sum(((test_ct(test_ct < boundry) - mean(test_ct(test_ct < boundry))).^2))));
%         lin_regb_r2_low = 1 - ((sum((lin_regb_low_pred - test_ct(test_ct < boundry)).^2))/(sum(((test_ct(test_ct < boundry) - mean(test_ct(test_ct < boundry))).^2))));
%         ann_rega_r2_high = 1 - ((sum((ann_rega_high_pred - test_ct(test_ct > boundry)).^2))/(sum(((test_ct(test_ct > boundry) - mean(test_ct(test_ct > boundry))).^2))));
%         ann_regb_r2_high = 1 - ((sum((ann_regb_high_pred - test_ct(test_ct > boundry)).^2))/(sum(((test_ct(test_ct > boundry) - mean(test_ct(test_ct > boundry))).^2))));
%         lin_rega_r2_high = 1 - ((sum((lin_rega_high_pred - test_ct(test_ct > boundry)).^2))/(sum(((test_ct(test_ct > boundry) - mean(test_ct(test_ct > boundry))).^2))));
%         lin_regb_r2_high = 1 - ((sum((lin_regb_high_pred - test_ct(test_ct > boundry)).^2))/(sum(((test_ct(test_ct > boundry) - mean(test_ct(test_ct > boundry))).^2))));
% 
%         ann_r2_avg = (ann_rega_r2_low + ann_regb_r2_low + ann_rega_r2_high + ann_rega_r2_high)/4;
%         lin_r2_avg = (lin_rega_r2_low + lin_regb_r2_low + lin_rega_r2_high + lin_rega_r2_high)/4;
% 
% 
%         results = [results; boundry, lin_rega_low_gof.rmse, lin_rega_high_gof.rmse, lin_regb_low_gof.rmse, lin_regb_high_gof.rmse, ann_rega_low_rmse, ann_rega_high_rmse, ann_regb_low_rmse, ann_regb_high_rmse];
% 
%         ann_rmse_sum = ann_rega_low_rmse + ann_regb_low_rmse + ann_rega_high_rmse + ann_regb_high_rmse;
%         lin_rmse_sum = lin_rega_low_gof.rmse + lin_rega_high_gof.rmse + lin_regb_low_gof.rmse + lin_regb_high_gof.rmse;
%         if ann_rmse_sum < best_ann_rmse
%             best_ann_rmse = ann_rmse_sum;
%             best_ann = [ann_rega_low, ann_regb_low, ann_rega_high, ann_regb_high];
%             best_ann_boundry = boundry;
%         end
%         if lin_rmse_sum < best_lin_rmse
%             best_lin_rmse = lin_rmse_sum;
%             best_lin = {lin_rega_low, lin_regb_low, lin_rega_high, lin_regb_high};
%             best_lin_boundry = boundry;
%         end
%         if ann_r2_avg > best_ann_r2_avg
%             best_ann_r2_avg = ann_r2_avg;
%             best_r2_ann = [ann_rega_low, ann_regb_low, ann_rega_high, ann_regb_high];
%             best_r2_ann_boundry = boundry;
%         end
%         if lin_r2_avg > best_lin_r2_avg
%             best_lin_r2_avg = lin_r2_avg;
%             best_r2_lin = {lin_rega_low, lin_regb_low, lin_rega_high, lin_regb_high};
%             best_r2_lin_boundry = boundry;
%         end
% 
%     end
% 
% 
% end
% 
% figure()
% hold on
% plot(results(:,1), results(:,2), 'DisplayName', "ELT A - Linear Low Calib.")
% plot(results(:,1), results(:,3), 'DisplayName', "ELT A - Linear High Calib.")
% plot(results(:,1), results(:,4), 'DisplayName', "ELT B - Linear Low Calib.")
% plot(results(:,1), results(:,5), 'DisplayName', "ELT B - Linear High Calib.")
% legend()
% xlabel("Low/High Calib. Boundry [ppm]")
% ylabel("RMSE Calib. [ppm]")
% title("Linear Performance vs. Boundry")
% figure()
% hold on
% 
% plot(results(:,1), results(:,6), 'DisplayName', "ELT A - Network Low Calib.")
% plot(results(:,1), results(:,7), 'DisplayName', "ELT A - Network High Calib.")
% plot(results(:,1), results(:,8), 'DisplayName', "ELT B - Network Low Calib.")
% plot(results(:,1), results(:,9), 'DisplayName', "ELT B - Network High Calib.")
% legend()
% xlabel("Low/High Calib. Boundry [ppm]")
% ylabel("RMSE Calib. [ppm]")
% title("Network Performance vs. Boundry")
% 


%% ANN Regression
ann_rega = feedforwardnet([16, 16]);
ann_regb = feedforwardnet([16, 16]);
ann_rega = train(ann_rega, train_ca', train_ct');
ann_regb = train(ann_regb, train_cb', train_ct');

%% Score Regressions
ann_rega_pred = ann_rega(test_ca')';
ann_regb_pred = ann_regb(test_cb')';
lin_rega_pred = predict(lin_rega, test_ca);
lin_regb_pred = predict(lin_regb, test_cb);

ann_rega_r2 = 1 - ((sum((ann_rega_pred - test_ct).^2))/(sum(((test_ct - mean(test_ct)).^2))));
ann_regb_r2 = 1 - ((sum((ann_regb_pred - test_ct).^2))/(sum(((test_ct - mean(test_ct)).^2))));
lin_rega_r2 = 1 - ((sum((lin_rega_pred - test_ct).^2))/(sum(((test_ct - mean(test_ct)).^2))));
lin_regb_r2 = 1 - ((sum((lin_regb_pred - test_ct).^2))/(sum(((test_ct - mean(test_ct)).^2))));

lin_rega_rmse = lin_rega.RMSE;
lin_regb_rmse = lin_regb.RMSE;
ann_rega_rmse = sqrt(mean((ann_rega_pred - test_ct).^2));
ann_regb_rmse = sqrt(mean((ann_regb_pred - test_ct).^2));
%% Examine Regression Residuals
ann_rega_resid = ann_rega_pred - test_ct;
ann_regb_resid = ann_regb_pred - test_ct;
lin_rega_resid = lin_rega_pred - test_ct;
lin_regb_resid = lin_regb_pred - test_ct;

fig4 = figure();
hold on
subplot(2, 2, 1)
plot(test_ct, ann_rega_resid, 'o')
title("ANN ELT A")
ylabel("CO_2 Residuals [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend("RMSE:    "+ann_rega_rmse+"\newlineR^2:    "+ann_rega_r2)
grid on
subplot(2, 2, 2)
plot(test_ct, ann_regb_resid, 'o')
title("ANN ELT B")
ylabel("CO_2 Residuals [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend("RMSE:    "+ann_regb_rmse+"\newlineR^2:    "+ann_regb_r2)
grid on
subplot(2, 2, 3)
plot(test_ct, lin_rega_resid, 'o')
title("Linear ELT A")
ylabel("CO_2 Residuals [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend("RMSE:    "+lin_rega_rmse+"\newlineR^2:    "+lin_rega_r2)
grid on
subplot(2, 2, 4)
plot(test_ct, lin_regb_resid, 'o')
title("Linear ELT B")
ylabel("CO_2 Residuals [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend("RMSE:    "+lin_regb_rmse+"\newlineR^2:    "+lin_regb_r2)
grid on
sgtitle("Comparison of Linear and ANN Regression Residuals")

%% Examine Regression Performance

fig5 = figure();
subplot(2, 2, 1)
hold on
plot(test_ct, ann_rega_pred, 'o')
plot(test_ct, test_ct, '-')
title("ANN ELT A")
ylabel("ANN Response CO_2 [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend(["RMSE:    "+ann_rega_rmse+"\newlineR^2:    "+ann_rega_r2,'1:1 Fit'])
grid on
subplot(2, 2, 2)
hold on
plot(test_ct, ann_regb_pred, 'o')
plot(test_ct, test_ct, '-')
title("ANN ELT B")
ylabel("ANN Response CO_2 [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend(["RMSE:    "+ann_regb_rmse+"\newlineR^2:    "+ann_regb_r2,'1:1 Fit'])
grid on
subplot(2, 2, 3)
hold on
plot(test_ct, lin_rega_pred, 'o')
plot(test_ct, test_ct, '-')
title("Linear ELT A")
ylabel("Linear Response CO_2 [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend(["RMSE:    "+lin_rega_rmse+"\newlineR^2:    "+lin_rega_r2,'1:1 Fit'])
grid on
subplot(2, 2, 4)
hold on
plot(test_ct, lin_regb_pred, 'o')
plot(test_ct, test_ct, '-')
title("Linear ELT B")
ylabel("Linear Response CO_2 [ppm]")
xlabel("Picarro CO_2 [ppm]")
legend(["RMSE:    "+lin_regb_rmse+"\newlineR^2:    "+lin_regb_r2,'1:1 Fit'])
grid on
sgtitle("Comparison of Linear and ANN Regression Performance")





fontsize(fig5, 30, 'points')
fontname(fig5, 'Times New Roman')



fig = figure();
hold on

plot(test_ct, test_ct, 'r--', 'LineWidth', 2);
plot(test_ct, ann_rega_pred, 'b.', 'MarkerSize', 20);

txt = "RMSE: " + ann_rega_rmse + " ppm\newlineR^2: " + ann_rega_r2;

text(min(xlim)+5, max(ylim)-40,  txt,'Interpreter','tex');
xlabel("ELT A CO_2 [ppm]",'Interpreter','tex');
ylabel("LICOR CO_2 [ppm]", 'Interpreter','tex');
title('Feed Forward Network Calibrations','Interpreter','tex');
legend(["1:1 Fit","Fitted CO_2 Dataset"], 'Interpreter', 'tex');
fontsize(fig,50, 'points')
fontname(fig, 'Times New Roman')
%%

fig = figure();
subplot(1, 2, 1)
hold on

plot(test_ct, test_ct, '--', 'LineWidth', 2);
plot(test_ct, ann_rega_pred, '.', 'MarkerSize', 20);

txt = "RMSE: " + round(ann_rega_rmse,3) + " ppm\newlineR^2: " + round(ann_rega_r2,3);

text(min(xlim)+5, max(ylim)-40,  txt,'Interpreter','tex');
ylabel("ELT A CO_2 [ppm]",'Interpreter','tex');
xlabel("Picarro CRDS CO_2 [ppm]", 'Interpreter','tex');
title('Performance','Interpreter','tex');
legend(["1:1 Fit","Fitted CO_2 Dataset"], 'Interpreter', 'tex');



daq = retime(daq, 'minutely');

times = daq.T < datetime(2024, 5, 4, 23,59,59,0) & daq.T > datetime(2024, 4, 28, 23,59,59,0);
daq = daq(times, :);
times = picarro.T < datetime(2024, 5, 4, 23,59,59,0) & picarro.T > datetime(2024, 4, 28, 23,59,59,0);
picarro = picarro(times, :);

subplot(1, 2, 2)
hold on
plot(picarro.T, picarro.CO2_sync, '.','LineWidth', 2.5, 'MarkerSize', 2)
plot(daq.T, daq.CA,'s', 'LineWidth', 2.5, 'MarkerSize', 2)
plot(daq.T, ann_rega([daq.CA, daq.TA, daq.HA]')','d', 'LineWidth', 2.5, 'MarkerSize', 2)
legend(["Picarro CRDS","NDIR Sensor", "Network Response (NDIR Sensor)"])
xlabel("Datetime")
ylabel("CO_2 [ppm]")
title("Performance")
grid on


sgtitle("Ambient Sensor Calibration")


fontsize(fig,50, 'points')
fontname(fig, 'Times New Roman')

%%
%save("../flux_test/calib", 'lin_rega', 'lin_regb', 'ann_rega', 'ann_regb')