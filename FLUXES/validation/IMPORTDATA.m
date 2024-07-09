function [daq,licor] = IMPORTDATA(path, dataset)
%IMPORTDATA Summary of this function goes here
    %   Detailed explanation goes here


    % import licor and daq dataset
    licor = IMPORTLICORFILE(path+dataset+'/licor.data');
    daq = IMPORTDAQFILE(path+dataset+'/daq.txt');

    % convert Q from ccm to lpm
    daq.Q = daq.Q/1000;
end

