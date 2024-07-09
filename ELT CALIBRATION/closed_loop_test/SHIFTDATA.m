function data = SHIFTDATA(data, lag)
%SHIFTDATA Shifts timetable based on lag amount.
%   This function shifts a dataset by a specified lag based on cross
%   corelation auto adjustment in the closed_loop.m function.
    if lag > 0
        data = [nan(lag, 1); data(1:end-lag)];
    elseif lag < 0
        data = [data(-lag+1:end); nan(-lag, 1)];
    end
    data = rmmissing(data);
end

