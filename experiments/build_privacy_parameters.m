clear all;
close all;
clc;

data = readtable("pubmed_probabilistic_data.csv");
unique_counts = readtable("pubmed_unique_counts.csv");
[unique_epsilon_vals, props] = reading_distribution(data);

all_rows = [];

for i = 2:length(props) % loop was started from 0 as the first column contained the user number, which is irrelevent. If the user number was not converted in the reading 
    % distribution function for a custom dataset, this value might need to
    % be modified. 
    [a_values, b_values, alpha_values, p_values, q_values] = privacy_parameters(unique_epsilon_vals{i}', props{i}', unique_counts{i-1,:}(2));
    
    %disp(alpha_values)
    % Ensure row vectors
    a_values = a_values(:)';
    b_values = b_values(:)';
    p_values = p_values(:)';
    q_values = q_values(:)';
    row = [a_values, b_values, alpha_values, p_values, q_values];

    % Append to master matrix
    all_rows = [all_rows; row];
end

% Define headers
headers = {'a1', 'a2', 'a3', 'b1', 'b2', 'b3', ...
           'alpha1',...
           'p1', 'p2', 'p3', 'q1', 'q2', 'q3'};

% Write to CSV
csv_filename = 'privacy_parameters_output.csv';
% I'm keping the outputs to download folder 
fid = fopen(csv_filename, 'w');
fprintf(fid, '%s,', headers{1:end-1});
fprintf(fid, '%s\n', headers{end});
fclose(fid);

% Append data
dlmwrite(csv_filename, all_rows, '-append');
disp("DONE CONVERTING TO CSV")
