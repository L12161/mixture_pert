%% Plot setting
clc; clear; close all;
set(0,'defaulttextinterpreter','latex'); 
set(0,'defaultlinelinewidth',1.3); 
set(0,'DefaultLineMarkerSize',6); 
set(0,'DefaultTextFontSize', 14); 
set(0,'DefaultAxesFontSize',14); 

figure;
h = plot([1 2],rand(6,2));
c = get(h,'Color');
Color = cell2mat(c);
Color = Color([1 2 5 3 4],:);
Color(3,:) = [0.13,0.55,0.13];
close all;

%% USER INPUT SECTION
% Ask user for epsilon values
fprintf('Enter the 3 epsilon values (space-separated): ');
epsilon_input = input('', 's');
epsilon_vals = str2num(epsilon_input);

if length(epsilon_vals) ~= 3
    error('Please provide exactly 3 epsilon values');
end

% Ask for epsilon values to skip
fprintf('Enter epsilon values to skip (space-separated, or press Enter for none): ');
skip_input = input('', 's');
if isempty(skip_input)
    epsilon_to_skip = [];
else
    epsilon_to_skip = str2num(skip_input);
end

% Ask for CSV file path and domain size
epsilon_csv = input('Enter path to epsilon CSV file (user,epsilon format): ', 's');
N_loc = input('Enter domain size: ');
selected_feature = 'Feature';

%% LOAD CSV FILE
fprintf('Loading CSV file...\n');

% Load epsilon values CSV (user, epsilon format)
epsilon_data = readtable(epsilon_csv);
feature_epsilons = epsilon_data.epsilon; % Read epsilon column
N_user = length(feature_epsilons);

fprintf('Loaded epsilon data: %d users\n', N_user);
fprintf('Domain size: %d\n', N_loc);

%% DERIVE EPSILON DISTRIBUTION
% Find unique epsilon values and their counts
[unique_eps, ~, idx] = unique(feature_epsilons);
epsilon_counts = accumarray(idx, 1);

% Calculate probabilities
Prob = epsilon_counts / N_user;

fprintf('Epsilon distribution:\n');
for i = 1:length(unique_eps)
    fprintf('  %.3f: %.1f%% (%d users)\n', unique_eps(i), Prob(i)*100, epsilon_counts(i));
end

% Map user epsilon values to levels (1, 2, 3)
N_lev = length(unique_eps);
if N_lev ~= 3
    error('Expected exactly 3 unique epsilon values, found %d', N_lev);
end

% Create mapping from epsilon values to levels
W = unique_eps' ./ epsilon_vals(1); % Normalize by first epsilon value
alpha = N_loc * Prob; 

%% CREATE W_LIST (privacy level for each location)
W_list = [];
for j = 1:N_lev
    W_list = [W_list; j*ones(round(alpha(j)),1)];
end     

% Adjust if rounding caused length mismatch
while length(W_list) < N_loc
    W_list = [W_list; N_lev]; % Add highest privacy level
end
W_list = W_list(1:N_loc); % Trim to exact size

W_list = W_list(randperm(N_loc));   
W_tab = tabulate(W_list);

W_index = cell(N_lev,1);
for j = 1:N_lev
    W_index{j} = find(W_list==j);
end

%% GENERATE SYNTHETIC DATA
data = generate_gaussian(N_user,N_loc);

% ===============  DATASET EVALUATION =======================
unique_vals = unique(data);              
counts = histc(data, unique_vals);       
mean_count = mean(counts);
std_count = std(counts);
fprintf ('Total number of unique counts: %.4f\n', length(unique_vals));
fprintf('Mean of counts: %.4f\n', mean_count);
fprintf('Standard deviation of counts: %.4f\n', std_count);
% ===============  DATASET EVALUATION =======================

%% OPTIMIZATION FUNCTIONS
fun0 = @(x) alpha'*( (x(N_lev+1:end) - x(N_lev+1:end).^2) ./ ( (x(1:N_lev)-x(N_lev+1:end)).^2 ) )...
    + max( (1-x(1:N_lev)-x(N_lev+1:end)) ./ (x(1:N_lev)-x(N_lev+1:end)) );

fun1 = @(x) alpha'*(exp(x)./((exp(x)-1).^2)); % symmetric
fun2 = @(b) alpha'*((b-b.^2)./((0.5-b).^2)) + 1; % a = 0.5

%% MAIN COMPUTATION LOOP
myEpsilon = epsilon_vals(1) * [0.25:0.25:6]; % Scale based on first epsilon value
N_epsilon = length(myEpsilon);

% Remove epsilon = 1 from the list
epsilon_to_skip = 1.0;
skip_indices = find(abs(myEpsilon - epsilon_to_skip) < 0.01); % Use small tolerance for floating point comparison
myEpsilon(skip_indices) = [];
N_epsilon = length(myEpsilon);

fprintf('Epsilon values to plot: [%s] (skipping %.2f)\n', num2str(myEpsilon), epsilon_to_skip);

alpha_values = [];
a_values = [];
b_values = [];
p_values = [];
q_values = [];

for i = 1:N_epsilon
    epsilon = myEpsilon(i);
    temp = W*epsilon;
    
    % Convert MSE to L1 distance using sqrt approximation
    mse_val_RAPPO = N_loc*exp(epsilon/2)./((exp(epsilon/2)-1).^2);
    result(1,i) = sqrt(mse_val_RAPPO * N_loc);
    a = ones(N_lev,1)*exp(epsilon/2)/(exp(epsilon/2)+1);
    b = 1-a;
    [L1_dist(1,i), ~] = actual_L1(a,b,data,N_loc,W_list);
    
    % result 2 deals with OUE - Convert MSE to L1 distance
    mse_val_OUE = N_loc*4*exp(epsilon)./((exp(epsilon)-1).^2)+1;
    result(2,i) = sqrt(mse_val_OUE * N_loc);
    a = ones(N_lev,1)*0.5; 
    b = ones(N_lev,1)/(exp(epsilon)+1);
    [L1_dist(2,i), ~] = actual_L1(a,b,data,N_loc,W_list);
    
    % result 3 deals with GRR - Convert MSE to L1 distance
    mse_val_GRR = N_loc * (N_loc - 2 + exp(epsilon)) / (exp(epsilon) - 1)^2;
    result(3,i) = sqrt(mse_val_GRR * N_loc);
    a = ones(N_lev,1) * exp(epsilon) / (exp(epsilon) + N_loc - 1);
    b = ones(N_lev,1) * 1 / (exp(epsilon) + N_loc - 1);
    [L1_dist(3,i), ~] = actual_L1(a,b,data,N_loc,W_list);

    % Convert optimization results from MSE to L1 distance
    [X, mse_result] = min_opt0(temp,fun0);     
    result_min(1,i) = sqrt(mse_result * N_loc);  
    Xmin0(:,i) = X;  
    a = X(1:N_lev);
    b = X((N_lev+1):end);
    [L1_dist_min(1,i), ~] = actual_L1(a,b,data,N_loc,W_list);  

    % fun3 with mixture model
    fun3 = @(x) x(4*N_lev+1)*sum((x(2*N_lev+1:3*N_lev).*(1-x(2*N_lev+1:3*N_lev)))./(x(2*N_lev+1:3*N_lev)-x(3*N_lev+1:4*N_lev)).^2) ...
        + sum(alpha(:) .* (1-x(4*N_lev+1)) .* x(N_lev+1:2*N_lev).*(1-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)).^2) ...
        + max((1-x(4*N_lev+1)) .* (1-x(1:N_lev)-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)));
    
    [X, mse_result_MIX] = min_opt3(temp,fun3, (numel(unique(data))));     
    result_min(4,i) = sqrt(mse_result_MIX * N_loc);  
    a = X(1:N_lev);
    b = X(N_lev+1:2*N_lev);
    p = X(2*N_lev+1:3*N_lev);
    q = X(3*N_lev+1:4*N_lev);
    alpha_param = X(4*N_lev+1);
    
    a_values = [a_values;a];
    b_values = [b_values;b];
    alpha_values = [alpha_values;alpha_param];
    p_values = [p_values; p];
    q_values = [q_values; q];
    
    % Use unbiased mixture estimator
    %[L1_dist_min(4,i), ~] = actual_L1_mixture(a,b,p,q,alpha_param,data,N_loc,W_list);
    
    [X, mse_result] = min_opt1(temp,fun1); 
    result_min(2,i) = sqrt(mse_result * N_loc);  
    Xmin1(:,i) = X;
    a = exp(X)./(1+exp(X)); 
    b = 1./(1+exp(X));
    [L1_dist_min(2,i), ~] = actual_L1(a,b,data,N_loc,W_list);
    
    [X, mse_result] = min_opt2(temp,fun2); 
    result_min(3,i) = sqrt(mse_result * N_loc);  
    Xmin2(:,i) = X;
    a = 0.5*ones(N_lev,1); 
    b = X;
    [L1_dist_min(3,i), ~] = actual_L1(a,b,data,N_loc,W_list);
end

%% ANOMALY DETECTION AND SMOOTHING
% Apply smoothing to all result arrays
fprintf('\nChecking for anomalies in results...\n');
result_smooth = result;
result_min_smooth = result_min;
L1_dist_smooth = L1_dist;
L1_dist_min_smooth = L1_dist_min;

for method = 1:size(result, 1)
    result_smooth(method, :) = fix_anomalies(result(method, :));
end

for method = 1:size(result_min, 1)
    result_min_smooth(method, :) = fix_anomalies(result_min(method, :));
end

for method = 1:size(L1_dist, 1)
    L1_dist_smooth(method, :) = fix_anomalies(L1_dist(method, :));
end

for method = 1:size(L1_dist_min, 1)
    L1_dist_min_smooth(method, :) = fix_anomalies(L1_dist_min(method, :));
end

%% PLOTTING
figure('Position',[100,100,1000,1000]);
ah1 = TightPlots(1, 1, 800,[10 7],[80,80],[50,50,60,40],[70,10],'pixels');
axes(ah1(1));

xlim_min = min(myEpsilon);  xlim_max = max(myEpsilon); 

% Calculate adaptive y-limits based on actual data
all_y_data = [result_smooth(1,:), result_smooth(2,:), result_smooth(3,:), result_min_smooth(4,:)];
min_y = min(all_y_data);
max_y = max(all_y_data);

% Add some padding (10% on each side)
y_range = max_y - min_y;
ylim_min = max(1, min_y - 0.1*y_range); % Don't go below 1 for log scale
ylim_max = max_y + 0.1*y_range;

fprintf('Y-axis limits: %.2f to %.2f\n', ylim_min, ylim_max);

% Draw background first (for x < median)
median_eps = median(myEpsilon);
fill([xlim_min median_eps median_eps xlim_min], [ylim_min ylim_min ylim_max ylim_max], ...
     [0.8 0.92 1], 'EdgeColor', 'none','FaceAlpha', 0.5); 
hold on;  

% Plot lines for L1 Distaance
h1=plot(myEpsilon, result(1,:), '-.o', 'Color', Color(1,:));
h2=plot(myEpsilon, result(2,:), '-.s', 'Color', Color(2,:));
h3=plot(myEpsilon, result(3,:), '-.^', 'Color', Color(3,:));
h4=plot(myEpsilon, result_min(4,:), '-.p', 'Color', Color(5,:));


% Plot lines for MSE Distaance
h1=plot(myEpsilon, mse_val_RAPPO, '-.o', 'Color', Color(1,:));
h2=plot(myEpsilon, mse_val_OUE, '-.s', 'Color', Color(2,:));
h3=plot(myEpsilon, mse_val_GRR, '-.^', 'Color', Color(3,:));
h4=plot(myEpsilon, mse_result_MIX, '-.p', 'Color', Color(5,:));

% Axes and plot settings
xlim([xlim_min, xlim_max]);
ylim([ylim_min, ylim_max]);
set(gca, 'YScale', 'log');

% Set adaptive y-ticks based on the data range
if ylim_max/ylim_min > 100
    % Wide range - use powers of 10
    yticks([1 10 100 1000 10000 100000]);           
else
    % Narrow range - use more granular ticks
    log_min = floor(log10(ylim_min));
    log_max = ceil(log10(ylim_max));
    ytick_values = 10.^(log_min:log_max);
    yticks(ytick_values);
end 
xlabel('$\epsilon$', 'FontSize', 35);
ylabel('L1 Distance', 'FontSize', 20);
legend([h1 h2 h3 h4], {'RAPPOR','OUE','GRR','Mixture'});
grid on;

% Add annotations with feature info
text('Parent',ah1(1),'string',{sprintf('Domain Size: %d', N_loc)},...
    'Position',[xlim_min + 0.6*(xlim_max-xlim_min), 18.4, 0]);

annotation('textbox',...
    [0.386 0.839909089283511 0.10950000231266 0.0823863652619449],...
    'String',{'High','Privacy'}, 'EdgeColor','none');
annotation('textbox',...
    [0.549 0.217045452919874 0.106500002241135 0.0823863652619449],...
    'String',{'Low','Privacy'}, 'EdgeColor','none');

% Save figure with feature name
filename = sprintf('L1_Distance_Analysis.png');
print(gcf, filename, '-dpng', '-r600');
fprintf('Saved plot as: %s\n', filename);

%% FUNCTIONS
function [L1_dist, RE] = actual_L1(a_list, b_list, data, N_loc, W_list)
repeat = 10;
N_user = length(data);

tab_true = tabulate(data);
count_true = tab_true(:,2);
L1_dist = zeros(1,repeat);
RE = zeros(N_loc,repeat);

parfor i = 1:repeat
    count_est = est(W_list,N_loc,N_user,a_list,b_list,count_true);
    L1_dist(i) = sum(abs(count_est-count_true));  
    RE(:,i) = abs(count_est-count_true)./count_true;
end

L1_dist = mean(L1_dist);
RE = mean(RE,2);
end

function est_count = est(W_list,N_loc,N_user,a_list,b_list,count_true)
est_count = zeros(N_loc,1);

parfor k = 1:N_loc
    n1 = count_true(k);
    n2 = N_user - n1;
    a = a_list(W_list(k)); 
    b = b_list(W_list(k)); 
    y1 = randsrc(n1,1,[[0 1]; [1-a a]]);
    y2 = randsrc(n2,1,[[0 1]; [1-b b]]);
    raw_count = sum([y1;y2]);

    est_count(k) = (raw_count-N_user*b)/(a-b);
end
end

% Function for unbiased mixture estimator
function [L1_dist, RE] = actual_L1_mixture(a_list, b_list, p_list, q_list, alpha_param, data, N_loc, W_list)
repeat = 10;
N_user = length(data);

tab_true = tabulate(data);
count_true = tab_true(:,2);

L1_dist = zeros(1,repeat);
RE = zeros(N_loc,repeat);

parfor i = 1:repeat
    count_est = est_mixture_simulated(W_list, N_loc, N_user, a_list, b_list, ...
                                      p_list, q_list, alpha_param, count_true);
    L1_dist(i) = sum(abs(count_est - count_true));
    RE(:,i) = abs(count_est - count_true) ./ count_true;
end

L1_dist = mean(L1_dist);
RE = mean(RE, 2);
end

function est_count = est_mixture_simulated(W_list, N_loc, N_user, a_list, b_list, ...
                                          p_list, q_list, alpha_param, count_true)
est_count = zeros(N_loc, 1);

parfor k = 1:N_loc
    n1 = count_true(k);
    n2 = N_user - n1;
    
    a = a_list(W_list(k));
    b = b_list(W_list(k));
    p = p_list(W_list(k));
    q = q_list(W_list(k));
    
    % Generate responses from mixture mechanism
    rand_vals_1 = rand(n1, 1);
    uses_grr_1 = rand_vals_1 < alpha_param;
    
    y1 = zeros(n1, 1);
    y1(uses_grr_1) = randsrc(sum(uses_grr_1), 1, [[0 1]; [1-p p]]);
    y1(~uses_grr_1) = randsrc(sum(~uses_grr_1), 1, [[0 1]; [1-a a]]);
    
    rand_vals_2 = rand(n2, 1);
    uses_grr_2 = rand_vals_2 < alpha_param;
    
    y2 = zeros(n2, 1);
    y2(uses_grr_2) = randsrc(sum(uses_grr_2), 1, [[0 1]; [1-q q]]);
    y2(~uses_grr_2) = randsrc(sum(~uses_grr_2), 1, [[0 1]; [1-b b]]);
    
    raw_count = sum([y1; y2]);
    
    numerator = raw_count - N_user * (alpha_param * q + (1 - alpha_param) * b);
    denominator = alpha_param * (p - q) + (1 - alpha_param) * (a - b);
    
    est_count(k) = numerator / denominator;
end
end

% Function to detect and fix anomalous points
function smoothed_data = fix_anomalies(data, threshold_factor)
    if nargin < 2
        threshold_factor = 3; % Default threshold
    end
    
    smoothed_data = data;
    n = length(data);
    
    if n < 3
        return; % Need at least 3 points
    end
    
    % Calculate differences between consecutive points
    diffs = abs(diff(data));
    
    % Calculate threshold based on median of differences
    median_diff = median(diffs);
    threshold = threshold_factor * median_diff;
    
    % Find anomalous points (skip first and last points)
    for i = 2:n-1
        % Check if current point creates large jumps on both sides
        left_jump = abs(data(i) - data(i-1));
        right_jump = abs(data(i+1) - data(i));
        
        % If both jumps are large, it's likely an anomaly
        if left_jump > threshold && right_jump > threshold
            % Additional check: see if removing this point creates smoother curve
            expected_value = (data(i-1) + data(i+1)) / 2;
            if abs(data(i) - expected_value) > threshold
                fprintf('Detected anomaly at point %d: %.2f -> %.2f\n', i, data(i), expected_value);
                smoothed_data(i) = expected_value;
            end
        end
    end
end