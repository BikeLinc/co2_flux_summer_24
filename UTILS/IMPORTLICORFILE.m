function licordata = IMPORTLICORFILE(filename)
%IMPORTLICORFILE Import data from a text file
%  DAQ = IMPORTLICORFILE(FILENAME) reads data from text file FILENAME for the
%  default selection.  Returns the data as a table.
%
%  DAQ = IMPORTLICORFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 22);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "DATE", "TIME", "H2O", "C", "CH4", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22"];
opts.SelectedVariableNames = ["DATE", "TIME", "C"];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "datetime", "datetime", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "DATE", "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, "TIME", "InputFormat", "HH:mm:ss");

% Import the data
licordata = readtable(filename, opts);

licordata.T = datetime( licordata.DATE + timeofday(licordata.TIME) , 'Format', 'default');
licordata.DATE = [];
licordata.TIME = [];
licordata = rmmissing(licordata);
licordata = table2timetable(licordata);

end
