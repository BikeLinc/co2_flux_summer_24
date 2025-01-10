function daq = IMPORTDAQFILE(filename, dataLines)
%IMPORTDAQFILE Imports a CSV-like file with no headers. Based on the number 
%of columns, variable names are assigned. The first column ("T") is always 
%a datetime variable, possibly in one of three formats. The remaining 
%columns are numeric.
%
%   DAQ = IMPORTDAQFILE(FILENAME) reads the file FILENAME for all rows 
%   (1, Inf). Returns a timetable.
%
%   DAQ = IMPORTDAQFILE(FILENAME, DATALINES) reads data for the specified 
%   row interval(s). 
%
%   Possible column sets:
%       7 columns : T, CA, TA, HA, CB, TB, HB
%       8 columns : T, Q, CA, TA, HA, CB, TB, HB
%      11 columns : T, MS, Q, CA, TA, HA, CB, TB, HB, T_AMB, T_CHMB

    %% Handle input arguments
    if nargin < 2
        dataLines = [1, Inf];
    end

    %% Step 1: Read file as a table with no headers
    opts = delimitedTextImportOptions();
    rawData = readtable(filename, opts);

    cols_non_empty = [];
    for colIdx = 1:width(rawData)
        if rawData{1, colIdx} ~= ""
            cols_non_empty = [cols_non_empty; colIdx];
        end
    end
    
    % Delete those columns from the original table
    rawData = rawData(:, cols_non_empty);

    %% Determine number of columns
    numCols = width(rawData);

    %% Step 2: Assign variable names based on number of columns
    switch numCols
        case 7
            varNames = ["T","CA","TA","HA","CB","TB","HB"];
        case 8
            varNames = ["T","Q","CA","TA","HA","CB","TB","HB"];
        case 11
            varNames = ["T","MS","Q","CA","TA","HA","CB","TB","HB","T_AMB","T_CHMB"];
        otherwise
            varNames = ["T","Q","CA","TA","HA","CB","TB","HB","T_AMB","T_CHMB"];
    end
    rawData.Properties.VariableNames = varNames;

    %% Step 3: Convert the first column ("T") to datetime
    % Attempt multiple known date/time formats in order
    dateFormats = { ...
        'yyyy-MM-dd''T''hh:mm:ss', ...
        'MM/dd/yy hh:mm:ss', ...
        'dd-MMM-yyyy HH:mm:ss', ...
        'MM/dd/yy hh:mm:ss a' ...
    };

    T_col = string(rawData.T);  % Convert to string for flexible parsing
    isParsed = false;
    for fmt = dateFormats
        try
            temp = datetime(T_col, 'InputFormat', fmt{1});
            if ~all(isnat(temp))   % If at least some entries were parsed
                rawData.T = datetime(T_col, 'InputFormat', fmt{1});
                isParsed = true;
                break;
            end
        catch
            % If parse fails for this format, try next
        end
    end

    % Get all variable names except 'T'
    vars = setdiff(rawData.Properties.VariableNames, 'T');
    
    % Loop through each variable to convert
    for v = vars
        varName = v{1};
        data = rawData.(varName);
        
        if iscell(data)
            rawData.(varName) = str2double(data);
        else
            rawData.(varName) = double(data);
        end
    end
    
    if ~isParsed
        error('Could not parse datetime column with known formats.');
    end

    %% Step 4: Convert other columns to numeric
    for n = 2:numCols  % columns 2..end
        rawData.(varNames{n}) = str2double(string(rawData.(varNames{n})));
    end

    %% Step 5: Convert to timetable
    daq = table2timetable(rawData, 'RowTimes', 'T');
end
