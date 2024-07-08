clc, clear, close all


%% Import Data

[file, location] = uigetfile("*.txt");

time = datetime(2024, 7, 5, 10, 52, 0);



%%
data = readtable(string(location) + string(file));

%%
data = renamevars(data,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8"],...
                ["MS","Q","CA","TA","HA","CB","TB","HB"]);

%% 

data.T = time + seconds(data.MS/1000);

%% writetable

writetable(data, string(location) +"_proc_"+ string(file))