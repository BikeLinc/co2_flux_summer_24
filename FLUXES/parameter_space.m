


% Data (For flux_measurement.m)
% temperature = daqTA;
% humidity = daqHA;
% ambientCO2 = daqCA;
% chamberCO2 = [daqCA; daqCB];

% Data (For calibration.m)
temperature = daq.TA;
humidity = daq.HA;
ambientCO2 = licor.C(licor.C>0);
chamberCO2 = 0;

% Fit normal distributions to the data
temp_mean = mean(temperature);
temp_std = std(temperature);
humid_mean = mean(humidity);
humid_std = std(humidity);
ambientCO2_mean = mean(ambientCO2);
ambientCO2_std = std(ambientCO2);
chamberCO2_mean = mean(chamberCO2);
chamberCO2_std = std(chamberCO2);

% Create a range of values for plotting
x_temp = linspace(min(temperature) - 2, max(temperature) + 2, 100);
x_humid = linspace(min(humidity) - 2, max(humidity) + 2, 100);
x_ambientCO2 = linspace(min(ambientCO2) - 10, max(ambientCO2) + 10, 100);
x_chamberCO2 = linspace(min(chamberCO2) - 10, max(chamberCO2) + 10, 100);

% Calculate the Gaussian distributions
pdf_temp = normpdf(x_temp, temp_mean, temp_std);
pdf_humid = normpdf(x_humid, humid_mean, humid_std);
pdf_ambientCO2 = normpdf(x_ambientCO2, ambientCO2_mean, ambientCO2_std);
pdf_chamberCO2 = normpdf(x_chamberCO2, chamberCO2_mean, chamberCO2_std);

% Create subplots
figure;

% --- Subplot 1: Temperature and Humidity ---
subplot(1, 2, 1);
hold on;
plot(x_temp, pdf_temp, 'b', 'LineWidth', 2); % Temperature plot in blue
area(x_temp, pdf_temp, 'FaceColor', 'b', 'FaceAlpha', 0.3); % Shading for Temperature
plot(x_humid, pdf_humid, 'g', 'LineWidth', 2); % Humidity plot in green
area(x_humid, pdf_humid, 'FaceColor', 'g', 'FaceAlpha', 0.3); % Shading for Humidity
xlabel('Value');
ylabel('Probability Density');
title('Temperature and Humidity Distribution');
legend('','T \circ C','', 'RH %');
hold off;

% --- Subplot 2: Ambient CO2 and Chamber CO2 ---
subplot(1, 2, 2);
hold on;
plot(x_ambientCO2, pdf_ambientCO2, 'r', 'LineWidth', 2); % Ambient CO2 plot in red
area(x_ambientCO2, pdf_ambientCO2, 'FaceColor', 'r', 'FaceAlpha', 0.3); % Shading for Ambient CO2
plot(x_chamberCO2, pdf_chamberCO2, 'm', 'LineWidth', 2); % Chamber CO2 plot in magenta
area(x_chamberCO2, pdf_chamberCO2, 'FaceColor', 'm', 'FaceAlpha', 0.3); % Shading for Chamber CO2
xlabel('Value');
ylabel('Probability Density');
title('Ambient and Chamber CO_2 Distribution');
legend('','Ambient CO_2','', 'Chamber CO_2');

sgtitle("Ambient Calibration Parameter Space")

hold off;

