% Define Parameters

Q_low       = 2.0833*10^-6;     % (m^3/s) -> .125 L/min
Q_high      = 6.6666*10^-5;     % (m^3/s) -> 4 L/min

uC_low      = 0.000415735;      % (mol/m^3) -> 10 ppm
uC_high     = 0.002078676;      % (mol/m^3) -> 50 ppm


F = (0.5:10)*10^-6;                   % 0.5 (umol/m^2/s)
uQ = 3.3333*10^-6;               % (m^3/s) -> 0.2 L/min
As = 0.063561;                   % (m^2)

n = 25;

Q_range = linspace(Q_low,Q_high,n);
uC_range = linspace(uC_low,uC_high,n);

[Q, uC] = meshgrid(Q_range, uC_range);

range = []


for idx = 1:length(F)
    F_idx =  F(idx);

    figure()
    uF = sqrt( ((F_idx./Q)*uQ).^2 +  2*((Q./As).*uC).^2 );

    
    [mz, idx] = min(uF(:))
    [xm, ym] = ind2sub(size(uF), idx);
    
    %axis tight;

    surf(Q*60000, uC, uF*10^6,'HandleVisibility','off')

    hold on;

    plot3(Q(xm, ym)*60000, uC(xm, ym), mz, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', "Optimal Flow Rate = " + Q(xm, ym)*60000 + " lpm");

    legend('Location', 'east');


    title("Flow Rate Optimization")
    xlabel("Flow Rate, Q [lpm]")
    ylabel("Enhancement Above Ambient [mol/m^3]")
    zlabel("Flux Rate Uncertainty [umol/m^2/s]")
    %%colorbar
    ax = gca;
    ax.ZAxis.Exponent = 0;
    
    minval = min(uF,[],'all');
    [row, col] = find(uF == minval);
    
    fprintf("Flux  : \t" + F_idx*10^6 + "\t umol/m^2/s")
    fprintf("Min-uF: \t" + minval*10^6 + "\t umol/m^2/s")
    fprintf("\t@ Q: \t" + Q(row,col)*60*1000 + "\t L/min")
    uC_ppm = (uC(row,col)*8.314*10^6*293.15)/(101325);
    fprintf("\t@ uC:  \t" + uC_ppm + "\t ppm")
    fprintf("Signal/Noise: \t" + F_idx/minval + "\t %%")
    
    range = [range; F_idx, minval, Q(row,col)*60*1000, uC_ppm, F_idx/minval];

    

end

range = array2table(range, "VariableNames", {'FLUX_MOL_M3_S', 'MIN_UNCERT_FLUX_MOL_M3_S', 'FLOW_L_MIN', 'UNCERT_CO2_PPM', 'SIGNAL_TO_NOISE'});

fig1 = figure()
plot(range.FLUX_MOL_M3_S*10^6, range.MIN_UNCERT_FLUX_MOL_M3_S*10^6)
title("Flux Rate vs. Uncertainty at Flux Rate")
xlabel("Flux Rate (\mu{}mol/m^2/s")
ylabel("Uncert. of Measuredment (\mu{}mol/m^2/s)")
grid on

fig2 = figure()
plot(range.FLUX_MOL_M3_S*10^6, range.FLOW_L_MIN)
title("Flux Rate vs. Uncertainty at Flux Rate")
xlabel("Flux Rate (\mu{}mol/m^2/s")
ylabel("Req. Flow Rate (L/min)")
grid on

fig3 = figure()
plot(range.FLUX_MOL_M3_S*10^6, range.SIGNAL_TO_NOISE)
title("Flux Rate vs. Uncertainty at Flux Rate")
xlabel("Flux Rate (\mu{}mol/m^2/s")
ylabel("Signal to Noise (%)")
grid on

fig4 = figure()
plot(range.FLOW_L_MIN, range.SIGNAL_TO_NOISE)
title("Flow Rate vs. Signal to Noise")
xlabel("Req. Flow Rate (L/min)")
ylabel("Signal to Noise (%)")
grid on



