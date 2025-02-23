function daqdata = IMPORTDAQ(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DAQDATA = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  DAQDATA = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  daqdata = importfile("S:\lascheer\CO2_FLUX\Test 4 April 2024\daq_data.txt", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 05-Apr-2024 09:37:35

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["T", "Q", "CA", "TA", "HA", "CB", "TB", "HB"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "T", "InputFormat", "yyyy-MM-dd'T'HH:mm:ss");

% Import the data
daqdata = readtable(filename, opts);

end