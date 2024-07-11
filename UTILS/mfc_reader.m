clc, clear, close all


mfc = serialport('COM1', 19200)

configureTerminator(mfc,"CR")

dataset = []
while 1
    writeline(mfc, "a");
    data = readline(mfc);
    data = strsplit(data);
    temp = double(data(2));
    pres = double(data(3));
    lpm = double(data(4));
    slpm = double(data(5));
   dataset = [dataset; temp pres lpm slpm]
   writematrix(dataset);
end



delete(mfc)