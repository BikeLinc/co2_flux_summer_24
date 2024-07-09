function data = SHIFTDATA(data, lag)
%SHIFTDATA Summary of this function goes here
%   Detailed explanation goes here
    if lag > 0
        data = [nan(lag, 1); data(1:end-lag)];
    elseif lag < 0
        data = [data(-lag+1:end); nan(-lag, 1)];
    end

end

