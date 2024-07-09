function [data, tss, Cchmb_tss] = CO2CHAMBERTSS(Cchmb_func, dt, Camb, F, As, Q, t0, V, uCchmb)
%MASSBALANCETSS Evaluates the CO2 Chamber mass balance function untill
%value reaches known steady state value within uncertainty range.

% calculate steady state value to compare with
Cchmb_tss = Camb + (F.*As)./Q;

% data will store computed values
data = zeros(1, 0);

% initializing loop variables
Cchmb = 0;
t = t0;

% run time to steady state calculation
while Cchmb <= abs(Cchmb_tss - uCchmb)
    Cchmb = Camb + (F.*As./Q).*(1-exp(-Q.*seconds(t-t0)./V));
    data = [data; [seconds(t), Cchmb]];
    t = t + dt;
end


headers = ["TIME", "CO2"];
data = array2table(data,'VariableNames', headers);

tss = t;

end

