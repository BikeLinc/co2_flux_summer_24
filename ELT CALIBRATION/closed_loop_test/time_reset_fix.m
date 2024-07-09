%% time_reset_fix.m
% Fixes issues where DAQ datasets have persistant time resetting to a
% constant value. This is because the microcontroller resets and the
% timestamp is not based on RTC value but is based on microcontroller
% millisecond value. (Retiring this method of timekeeping, incorrect)
%
% Lincoln Scheer
% 7/8/2024
%
%

clc, clear, close all

% Correct the timestamps
daq = IMPORTDAQFILE("data/7.5.2024/daq.csv");
daq = rmmissing(daq);



last_ts = daq.T(1);
for i = 2:height(daq)
    ts = daq.T(i);
    dts = seconds(ts - last_ts);
    if abs(dts) > 4
        disp(['Reset detected at index: ', num2str(i)]);
        for j = i:height(daq)
            daq.T(j) = last_ts + seconds(3);
            last_ts = daq.T(j);
        end
        break;
    end
    last_ts = ts;
end

writetimetable(daq, "data/7.5.2024/daq_TS_FIX.csv")


%%
figure()
hold on
plot(daq.T, daq.CA)
