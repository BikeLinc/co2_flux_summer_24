function daq = IMPORTDAQFILE(filename)
%IMPORTFILE Import data from a text file
%  DAQ = IMPORTDAQFILE(FILENAME) reads data from text file FILENAME for the
%  default selection.  Returns the data as a table.
%
%  DAQ = IMPORTDAQFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%

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
opts.SelectedVariableNames = ["T", "Q", "CA", "TA", "HA", "CB", "TB", "HB"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties06/16/24 11:18:26
%opts = setvaropts(opts, "T", "InputFormat", "MM/dd/yy HH:mm:ss");
opts = setvaropts(opts, "T", "InputFormat", 'yyyy-MM-dd''T''HH:mm:ss');

% Import the data
daq = readtable(filename, opts);
daq = table2timetable(daq);

end