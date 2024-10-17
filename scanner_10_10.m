%% Main Script
clc, clear, close all

% Define the directory to search
searchDir = pwd; % Current directory, change as needed

% Define the file extensions to search for
fileExtensions = {'.csv', '.mat', '.txt', '.TXT', '.data'}; % Add more extensions as needed

% Search for data files
dataFiles = searchForDataFiles(searchDir, fileExtensions);

% Display the found data files
disp('Found data files:');
disp(dataFiles);

% Search for empty files
emptyFiles = searchForEmptyFiles(searchDir);

% Display the found empty files
%disp('Found empty files:');
%disp(emptyFiles);

% Generate a report for empty files
generateEmptyFilesReport(emptyFiles);

% Find unused files
unusedFiles = findUnusedFiles(searchDir, dataFiles);

% Display the found unused files
disp('Found unused files:');
disp(unusedFiles);

% Generate a report for unused files
generateUnusedFilesReport(unusedFiles);

%% Function Definitions

function dataFiles = searchForDataFiles(searchDir, fileExtensions)
    % Initialize an empty cell array to store file paths
    dataFiles = {};

    % Recursively search for files
    files = dir(fullfile(searchDir, '**', '*'));
    for i = 1:length(files)
        % Skip directories
        if files(i).isdir
            continue;
        end

        % Get the file extension
        [~, ~, ext] = fileparts(files(i).name);

        % Check if the file extension matches any of the specified extensions
        if ismember(ext, fileExtensions)
            % Add the file path to the list
            dataFiles{end+1} = fullfile(files(i).folder, files(i).name);
        end
    end

    % Generate a report as a txt file called directory_listing
    fid = fopen('directory_listing.txt', 'w');
    for i = 1:length(dataFiles)
        % Print file type, name, and parent directory
        fprintf(fid, 'File Type: %s\n', fileExtensions{1});
        fprintf(fid, 'File Name: %s\n', dataFiles{i});
        fprintf(fid, 'Parent Directory: %s\n', fileparts(dataFiles{i}));
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function emptyFiles = searchForEmptyFiles(searchDir)
    % Initialize an empty cell array to store file paths
    emptyFiles = {};

    % Recursively search for files
    files = dir(fullfile(searchDir, '**', '*'));
    for i = 1:length(files)
        % Skip directories
        if files(i).isdir
            continue;
        end

        % Check if the file size is zero
        if files(i).bytes == 0
            % Add the file path to the list
            emptyFiles{end+1} = fullfile(files(i).folder, files(i).name);
        end
    end
end

function generateEmptyFilesReport(emptyFiles)
    % Generate a report for empty files
    fid = fopen('empty_files_report.txt', 'w');
    for i = 1:length(emptyFiles)
        % Print file name and parent directory
        fprintf(fid, 'File Name: %s\n', emptyFiles{i});
        fprintf(fid, 'Parent Directory: %s\n', fileparts(emptyFiles{i}));
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function unusedFiles = findUnusedFiles(searchDir, dataFiles)
    % Initialize an empty cell array to store unused file paths
    unusedFiles = {};

    % Recursively search for .m files
    scriptFiles = dir(fullfile(searchDir, '**', '*.m'));

    % Check each data file to see if it is referenced in any script file
    for i = 1:length(dataFiles)
        isUsed = false;
        for j = 1:length(scriptFiles)
            % Read the content of the script file
            scriptContent = fileread(fullfile(scriptFiles(j).folder, scriptFiles(j).name));
            % Check if the data file name is referenced in the script content
            if contains(scriptContent, dataFiles{i})
                isUsed = true;
                break;
            end
        end
        % If the data file is not used in any script, add it to the unused files list
        if ~isUsed
            unusedFiles{end+1} = dataFiles{i};
        end
    end
end

function generateUnusedFilesReport(unusedFiles)
    % Generate a report for unused files
    fid = fopen('unused_files_report.txt', 'w');
    for i = 1:length(unusedFiles)
        % Print file name and parent directory
        fprintf(fid, 'File Name: %s\n', unusedFiles{i});
        fprintf(fid, 'Parent Directory: %s\n', fileparts(unusedFiles{i}));
        fprintf(fid, '\n');
    end
    fclose(fid);
end