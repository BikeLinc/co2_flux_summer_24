function [offset] = XCORRANALYSIS(A,B)
%XCORRANALYSIS Summary of this function goes here
%   Detailed explanation goes here
    
[r, lags] = xcorr(A, B);

figure()
plot(lags, r)



end

