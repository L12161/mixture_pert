%% Privacy Budget vs Envelope Gap Plotter
% This script reads a CSV file and plots curves for different methods
% against mean privacy budget

clear all; close all; clc;

% Prompt user for CSV file
[filename, filepath] = uigetfile('*.csv', 'Select CSV file');

% Check if user cancelled
if isequal(filename, 0)
    disp('File selection cancelled.');
    return;
end

% Full file path
fullpath = fullfile(filepath, filename);

% Read the CSV file
data = readtable(fullpath);

% Extract columns
mean_eps = data.mean_eps;
grr_ue = data.("GRR_UE");
grr_olh = data.("GRR_OLH");
ue_olh = data.("UE_OLH");

% Create figure
figure('Position', [100, 100, 1000, 600]);

% Plot the three curves
hold on;
plot(mean_eps, grr_ue, 'o-', 'LineWidth', 4, 'MarkerSize', 6, 'DisplayName', 'GRR+UE', 'Color', [0 0.4470 0.7410]);
plot(mean_eps, grr_olh, '--', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'GRR+OLH', 'Color', [0.8500 0.3250 0.0980]);
plot(mean_eps, ue_olh, '--', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'UE+OLH', 'Color', [0.9290 0.6940 0.1250]);
hold off;

% Set labels and title
xlabel('Mean Privacy Budget', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Max Envelope Gap to Oracle', 'FontSize', 16, 'FontWeight', 'bold');
title('Oracle-Tightness Across Mean Privacy Budget', 'FontSize', 18, 'FontWeight', 'bold');

% Add legend
legend('Location', 'best', 'FontSize', 11);

% Add grid for better readability
grid on;
grid minor;

% Adjust layout
set(gca, 'FontSize', 10);

% Create inset axes for magnified view (x range 5 to 10)
% Position: [left, bottom, width, height] in normalized units
axes('Position', [0.55, 0.45, 0.35, 0.35]);
hold on;
plot(mean_eps, grr_ue, 'o-', 'LineWidth', 3, 'MarkerSize', 5, 'DisplayName', 'GRR+UE', 'Color', [0 0.4470 0.7410]);
plot(mean_eps, grr_olh, '--', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'GRR+OLH', 'Color', [0.8500 0.3250 0.0980]);
plot(mean_eps, ue_olh, '--', 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'UE+OLH', 'Color', [0.9290 0.6940 0.1250]);
hold off;

% Set x-axis range to 5-10 for magnified view
xlim([5, 9]);

xlabel('Mean Privacy Budget', 'FontSize', 14);
ylabel('Max Envelope Gap to Oracle', 'FontSize', 14);
title('Magnified View', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 8);