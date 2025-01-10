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
    % Define conversion factor
    conversion_factor = 44.01 * 1e3 * 120; % Converts μmol/m²/s to mg/m²/hour

    % Convert Units
    Q_lpm = cfg.cms_to_lpm(Q); % m³/s -> L/min
    uC_ppm = cfg.mol_to_ppm(uC); % mol/m³ -> ppm

    % Convert uF and F_idx
    uF_mg = uF * conversion_factor; % μmol/m²/s -> mg/m²/hour
    F_idx_mg = F_idx * conversion_factor; % μmol/m²/s -> mg/m²/hour

    % Find minimum uncertainty
    [mz, idx] = min(uF_mg(:));
    [xm, ym] = ind2sub(size(uF_mg), idx);

    % Plotting
    figure();
    surf(Q_lpm, uC_ppm, uF_mg, 'HandleVisibility', 'off');
    hold on;

    % Plot optimal point with 2 significant figures
    plot3(Q_lpm(xm, ym), uC_ppm(xm, ym), mz, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf("Optimal Flow Rate = %.2g L/min\nSignal/Noise = %.2g%%\n+/- Uncertainty = %.2g mg/m²/hour", ...
        Q_lpm(xm, ym), F_idx_mg/mz, mz));

    legend('Location', 'east');
    title(sprintf("Flow Rate Optimization for F = %.2g mg/m²/hour", F_idx_mg));
    xlabel('Flow Rate (L/min)');
    ylabel('Uncertainty in CO₂ (ppm)');
    zlabel('Uncertainty in Flux (mg/m²/hour)');
    grid on;

    % Print Results
    fprintf("Flux  : \t%.2f mg/m²/hour\n", F_idx_mg);
    fprintf("Min-uF: \t%.2f mg/m²/hour\n", mz);
    fprintf("\t@ Q: \t%.2f L/min\n", Q_lpm(xm, ym));
    fprintf("\t@ uC:  \t%.2f ppm\n", uC_ppm(xm, ym));
    fprintf("Signal/Noise: \t%.2f %%\n", F_idx_mg/mz);
end


% Main loop to process each F value
for idx = 1:length(F)
    F_idx = F(idx);
    uF = calculateUncertainty(F_idx, Q, uQ, As, uC);
    plotResults(Q, uC, uF, F_idx, cfg);
end