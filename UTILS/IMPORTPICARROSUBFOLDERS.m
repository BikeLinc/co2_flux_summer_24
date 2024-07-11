function data = importpicarrosubfolders(folder)
%IMPORTPICARROSUBFOLDERS Concatonates picarro .dat datasets in subfolders
%and saves as a .csv (can be big) and returns as a timetable object.

picFolder   = folder;
picFiles = dir(fullfile(picFolder, '**', '*.dat'));

data = [];
m = length(picFiles);
for idx = 1:m
    file = string(picFiles(idx).folder) + "/" + string(picFiles(idx).name);
    data = [data;importpicarrofile(file)];
    waitbar(idx/m) % Update progress bar
end

writetimetable(data, 'picarro.txt');

end