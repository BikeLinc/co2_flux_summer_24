%% Flow Rate Optimization for Flux Measurement
clc, clear, close all

% Add Path to Utility Functions
addpath(genpath('UTILS'));
cfg = config();

% Define Parameters (Normal Units)
% Flow Rate Range
Q_low = 0.125; % L/min
Q_high = 4; % L/min

% Uncertainty in CO2 Range
uC_low = 10; % ppm
uC_high = 50; % ppm

% Uncertainty in Flow Rate (% of Reading)
uQ = 0.33; % L/min

% Flux Range
F = (0.5:10);    % umol/m^2/s

% Chamber Surface Area
As = cfg.As;          % (m^2)

% Convert Units to SI
Q_low = cfg.lpm_to_cms(Q_low); % (m^3/s)
Q_high = cfg.lpm_to_cms(Q_high); % (m^3/s)
uQ = cfg.lpm_to_cms(uQ); % (m^3/s)
uC_low = cfg.ppm_to_mol(uC_low); % (mol/m^3)
uC_high = cfg.ppm_to_mol(uC_high); % (mol/m^3)
F = F * 10^-6; % (umol/m^2/s);

% Simulaiton Resolution (# of Points)
n = 25;

% Create Meshgrid
Q_range = linspace(Q_low, Q_high, n);
uC_range = linspace(uC_low, uC_high, n);
[Q, uC] = meshgrid(Q_range, uC_range);

% Function to calculate uncertainty
function uF = calculateUncertainty(F_idx, Q, uQ, As, uC)
    uF = sqrt(((F_idx ./ Q) * uQ).^2 + 2 * ((Q ./ As) .* uC).^2);
end

% Function to plot the results
function plotResults(Q, uC, uF, F_idx, cfg)

    % Convert Units (Q: m^3/s -> L/min, uC: mol/m^3 -> ppm, uF: umol/m^2/s -> mg/m^2/hour)
    Q = cfg.cms_to_lpm(Q);
    uC = cfg.mol_to_ppm(uC);
    % convert uF from mol to mg (MW = 44.01 g/mol) and from s to hour
    uF = uF * 44.01 * 10^3 * 60 * 60;
    F_idx = F_idx * 44.01 * 10^3 * 60 * 60;

    [mz, idx] = min(uF(:));
    [xm, ym] = ind2sub(size(uF), idx);
    

    % New Figure w/ Units from Conversions Above
    figure();
    surf(Q, uC, uF, 'HandleVisibility', 'off');
    hold on;


    % Repeat Above Dot But, Add Signal to Noise and +/- Uncert Value in Legend
    plot3(Q(xm, ym), uC(xm, ym), mz, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
        'DisplayName', "Optimal Flow Rate = " + Q(xm, ym) + " lpm" + newline + ...
        "Signal/Noise = " + F_idx/mz + "%" + newline + ...
        "+/- Uncertainty = " + mz + " mg/m^2/hour");



    legend('Location', 'east');
    title("Flow Rate Optimization for F = " + F_idx + " mg/m^2/hour");
    xlabel('Flow Rate (L/min)');
    ylabel('Uncertainty in CO2 (ppm)');
    zlabel('Uncertainty in Flux (mg/m^2/hour)');
    grid on;

    % Print Results
    fprintf("Flux  : \t" + F_idx + "\t mg/m^2/hour\n");
    fprintf("Min-uF: \t" + mz + "\t mg/m^2/hour\n");
    fprintf("\t@ Q: \t" + Q(xm, ym) + "\t L/min\n");
    fprintf("\t@ uC:  \t" + uC(xm, ym) + "\t ppm\n");
    fprintf("Signal/Noise: \t" + F_idx/mz + "\t %%\n");
end

% Main loop to process each F value
for idx = 1:length(F)
    F_idx = F(idx);
    uF = calculateUncertainty(F_idx, Q, uQ, As, uC);
    plotResults(Q, uC, uF, F_idx, cfg);
end