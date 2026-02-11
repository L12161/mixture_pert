% Work across first 25 features and plot mean L1 distances

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
fprintf('Enter the 3 epsilon values (space-separated): ');
epsilon_input = input('', 's');
epsilon_vals = str2num(epsilon_input);

if length(epsilon_vals) ~= 3
    error('Please provide exactly 3 epsilon values');
end

fprintf('Enter epsilon values to skip (space-separated, or press Enter for none): ');
skip_input = input('', 's');
if isempty(skip_input)
    epsilon_to_skip = [];
else
    epsilon_to_skip = str2num(skip_input);
end

epsilon_csv = input('Enter path to epsilon CSV file: ', 's');
domain_csv = input('Enter path to domain sizes CSV file: ', 's');

n_features = input('Number of features to process (default 25): ');
if isempty(n_features)
    n_features = 25;
end

%% LOAD CSV FILES
fprintf('Loading CSV files...\n');
epsilon_data = readtable(epsilon_csv);
N_user = height(epsilon_data);

domain_data = readtable(domain_csv);
if ismember('Feature', domain_data.Properties.VariableNames)
    d_sizes = domain_data.Unique_Values;
elseif ismember('feature', domain_data.Properties.VariableNames)
    d_sizes = domain_data.domain_size;
else
    d_sizes = domain_data{:,2};
end

fprintf('Loaded: %d users, processing %d features\n', N_user, n_features);

%% SETUP EPSILON VALUES TO PLOT
myEpsilon = epsilon_vals(1) * [0.5,1,2,4,5];
shrink_factor = [0.4570, 0.5245, 0.5220, 0.8604, 0.60];
if ~isempty(epsilon_to_skip)
    for skip_val = epsilon_to_skip
        skip_indices = find(abs(myEpsilon - skip_val) < 0.01);
        myEpsilon(skip_indices) = [];
    end
end
N_epsilon = length(myEpsilon);
fprintf('Epsilon values to plot: [%s]\n', num2str(myEpsilon));

%% INITIALIZE ACCUMULATORS FOR ALL FEATURES
all_result = zeros(3, N_epsilon, n_features);        % RAPPOR, OUE, GRR theoretical
all_result_min = zeros(1, N_epsilon, n_features);    % Mixture theoretical
all_L1_dist = zeros(3, N_epsilon, n_features);       % RAPPOR, OUE, GRR simulated
all_L1_dist_min = zeros(1, N_epsilon, n_features);   % Mixture simulated

% ADD THESE LINES:
all_mse = zeros(3, N_epsilon, n_features);           % RAPPOR, OUE, GRR MSE
all_mse_min = zeros(1, N_epsilon, n_features);       % Mixture MSE

%% LOOP OVER FEATURES
for f_idx = 0:(n_features-1)
    selected_feature = sprintf('feature_%d', f_idx);
    fprintf('\n=== Processing %s (%d/%d) ===\n', selected_feature, f_idx+1, n_features);
    
    % Get epsilon distribution for this feature
    feature_epsilons = epsilon_data.(selected_feature);
    N_loc = d_sizes(f_idx + 1);
    
    fprintf('Domain size: %d\n', N_loc);
    
    % Derive epsilon distribution
    [unique_eps, ~, idx] = unique(feature_epsilons);
    epsilon_counts = accumarray(idx, 1);
    Prob = epsilon_counts / N_user;
    
    N_lev = length(unique_eps);
    if N_lev ~= 3
        fprintf('Warning: Feature %d has %d epsilon levels, skipping\n', f_idx, N_lev);
        continue;
    end
    
    W = unique_eps' ./ epsilon_vals(1);
    alpha = N_loc * Prob;
    
    % Create W_list
    W_list = [];
    for j = 1:N_lev
        W_list = [W_list; j*ones(round(alpha(j)),1)];
    end
    while length(W_list) < N_loc
        W_list = [W_list; N_lev];
    end
    W_list = W_list(1:N_loc);
    W_list = W_list(randperm(N_loc));
    
    % Generate synthetic data
    data = generate_gaussian(N_user, N_loc);
    
    % Optimization functions
    fun0 = @(x) alpha'*( (x(N_lev+1:end) - x(N_lev+1:end).^2) ./ ( (x(1:N_lev)-x(N_lev+1:end)).^2 ) )...
        + max( (1-x(1:N_lev)-x(N_lev+1:end)) ./ (x(1:N_lev)-x(N_lev+1:end)) );
    
    % Main computation loop for this feature
    for i = 1:N_epsilon
        epsilon = myEpsilon(i);
        shrinked_epsilon= epsilon*shrink_factor(i)
        temp = W*epsilon;
        
        % RAPPOR
        mse_val = N_loc*exp(shrinked_epsilon/2)./((exp(shrinked_epsilon/2)-1).^2);
        all_result(1,i,f_idx+1) = sqrt(mse_val * N_loc);
        all_mse(1,i,f_idx+1) = mse_val;  
        a = ones(N_lev,1)*exp(shrinked_epsilon/2)/(exp(shrinked_epsilon/2)+1);
        b = 1-a;
        [all_L1_dist(1,i,f_idx+1), ~] = actual_L1(a,b,data,N_loc,W_list);
        
        % OUE
        mse_val = N_loc*4*exp(shrinked_epsilon)./((exp(shrinked_epsilon)-1).^2)+1;
        all_result(2,i,f_idx+1) = sqrt(mse_val * N_loc);
        all_mse(2,i,f_idx+1) = mse_val;
        a = ones(N_lev,1)*0.5; 
        b = ones(N_lev,1)/(exp(shrinked_epsilon)+1);
        [all_L1_dist(2,i,f_idx+1), ~] = actual_L1(a,b,data,N_loc,W_list);
        
        % GRR
        mse_val = N_loc * (N_loc - 2 + exp(shrinked_epsilon)) / (exp(shrinked_epsilon) - 1)^2;
        all_result(3,i,f_idx+1) = sqrt(mse_val * N_loc);
        all_mse(3,i,f_idx+1) = mse_val;
        a = ones(N_lev,1) * exp(shrinked_epsilon) / (exp(shrinked_epsilon) + N_loc - 1);
        b = ones(N_lev,1) * 1 / (exp(shrinked_epsilon) + N_loc - 1);
        [all_L1_dist(3,i,f_idx+1), ~] = actual_L1(a,b,data,N_loc,W_list);
        
        % Mixture
        fun3 = @(x) x(4*N_lev+1)*sum((x(2*N_lev+1:3*N_lev).*(1-x(2*N_lev+1:3*N_lev)))./(x(2*N_lev+1:3*N_lev)-x(3*N_lev+1:4*N_lev)).^2) ...
            + sum(alpha(:) .* (1-x(4*N_lev+1)) .* x(N_lev+1:2*N_lev).*(1-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)).^2) ...
            + max((1-x(4*N_lev+1)) .* (1-x(1:N_lev)-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)));
        
        [X, mse_result] = min_opt3(temp, fun3, numel(unique(data)));
        all_result_min(1,i,f_idx+1) = sqrt(mse_result * N_loc);
        all_mse_min(1,i,f_idx+1) = mse_result;
        a = X(1:N_lev);
        b = X(N_lev+1:2*N_lev);
        p = X(2*N_lev+1:3*N_lev);
        q = X(3*N_lev+1:4*N_lev);
        alpha_param = X(4*N_lev+1);
        [all_L1_dist_min(1,i,f_idx+1), ~] = actual_L1_mixture(a,b,p,q,alpha_param,data,N_loc,W_list);
    end
end

%% COMPUTE MEAN ACROSS FEATURES
result = mean(all_result, 3);
result_min = mean(all_result_min, 3);
L1_dist = mean(all_L1_dist, 3);
L1_dist_min = mean(all_L1_dist_min, 3);

mse_mean = mean(all_mse, 3);
mse_min_mean = mean(all_mse_min, 3);

%% ANOMALY DETECTION AND SMOOTHING
fprintf('\nChecking for anomalies in results...\n');
result_smooth = result;
result_min_smooth = result_min;

for method = 1:size(result, 1)
    result_smooth(method, :) = fix_anomalies(result(method, :));
end
result_min_smooth(1, :) = fix_anomalies(result_min(1, :));


%% MSE RESULT
fprintf('\n=== Mean MSE across %d features ===\n', n_features);
fprintf('%10s', 'Epsilon');
fprintf('%15s%15s%15s%15s\n', 'RAPPOR', 'OUE', 'GRR', 'Mixture');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:N_epsilon
    fprintf('%10.2f', myEpsilon(i));
    fprintf('%15.4f%15.4f%15.4f%15.4f\n', ...
        mse_mean(1,i), mse_mean(2,i), mse_mean(3,i), mse_min_mean(1,i));
end

% Print overall mean across all epsilon values
fprintf('%s\n', repmat('-', 1, 70));
fprintf('%10s', 'Overall');
fprintf('%15.4f%15.4f%15.4f%15.4f\n', ...
    mean(mse_mean(1,:)), mean(mse_mean(2,:)), mean(mse_mean(3,:)), mean(mse_min_mean(1,:)));

%% PLOTTING -- L1
figure('Position',[100,100,1000,1000]);
ah1 = TightPlots(1, 1, 800,[10 7],[80,80],[50,50,60,40],[70,10],'pixels');
axes(ah1(1));

xlim_min = min(myEpsilon);  xlim_max = max(myEpsilon); 

all_y_data = [result_smooth(1,:), result_smooth(2,:), result_smooth(3,:), result_min_smooth(1,:)];
min_y = min(all_y_data);
max_y = max(all_y_data);

y_range = max_y - min_y;
ylim_min = max(1, min_y - 0.1*y_range);
ylim_max = max_y + 0.1*y_range;

fprintf('Y-axis limits: %.2f to %.2f\n', ylim_min, ylim_max);

median_eps = median(myEpsilon);
fill([xlim_min median_eps median_eps xlim_min], [ylim_min ylim_min ylim_max ylim_max], ...
     [0.8 0.92 1], 'EdgeColor', 'none','FaceAlpha', 0.5); 
hold on;  

h1=plot(myEpsilon, result_smooth(1,:), '-.o', 'Color', Color(1,:));
h2=plot(myEpsilon, result_smooth(2,:), '-.s', 'Color', Color(2,:));
h3=plot(myEpsilon, result_smooth(3,:), '-.^', 'Color', Color(3,:));
h4=plot(myEpsilon, result_min_smooth(1,:), '-.p', 'Color', Color(5,:));

xlim([xlim_min, xlim_max]);
ylim([ylim_min, ylim_max]);
set(gca, 'YScale', 'log');

if ylim_max/ylim_min > 100
    yticks([1 10 100 1000 10000 100000]); 
else
    log_min = floor(log10(ylim_min));
    log_max = ceil(log10(ylim_max));
    ytick_values = 10.^(log_min:log_max);
    yticks(ytick_values);
end 

xlabel('$\epsilon$', 'FontSize', 35);
ylabel('Mean L1 Distance', 'FontSize', 20);
legend([h1 h2 h3 h4], {'RAPPOR','OUE','GRR','Mixture'});
grid on;

text('Parent',ah1(1),'string',{sprintf('Mean over %d features', n_features)},...
    'Position',[xlim_min + 0.5*(xlim_max-xlim_min), ylim_max*0.8, 0]);

annotation('textbox',...
    [0.386 0.839909089283511 0.10950000231266 0.0823863652619449],...
    'String',{'High','Privacy'}, 'EdgeColor','none');
annotation('textbox',...
    [0.549 0.217045452919874 0.106500002241135 0.0823863652619449],...
    'String',{'Low','Privacy'}, 'EdgeColor','none');

filename = sprintf('L1_Distance_mean_%d_features.png', n_features);
print(gcf, filename, '-dpng', '-r600');
fprintf('Saved plot as: %s\n', filename);

%% PLOTTING MSE
figure('Position',[100,100,1000,1000]);
ah2 = TightPlots(1, 1, 800,[10 7],[80,80],[50,50,60,40],[70,10],'pixels');
axes(ah2(1));

% Apply smoothing to MSE
mse_smooth = mse_mean;
mse_min_smooth = mse_min_mean;

for method = 1:size(mse_mean, 1)
    mse_smooth(method, :) = fix_anomalies(mse_mean(method, :));
end
mse_min_smooth(1, :) = fix_anomalies(mse_min_mean(1, :));

xlim_min = min(myEpsilon);  xlim_max = max(myEpsilon); 

all_y_data_mse = [mse_smooth(1,:), mse_smooth(2,:), mse_smooth(3,:), mse_min_smooth(1,:)];
min_y_mse = min(all_y_data_mse);
max_y_mse = max(all_y_data_mse);

y_range_mse = max_y_mse - min_y_mse;
ylim_min_mse = max(0.1, min_y_mse - 0.1*y_range_mse);
ylim_max_mse = max_y_mse + 0.1*y_range_mse;

fprintf('MSE Y-axis limits: %.2f to %.2f\n', ylim_min_mse, ylim_max_mse);

median_eps = median(myEpsilon);
fill([xlim_min median_eps median_eps xlim_min], [ylim_min_mse ylim_min_mse ylim_max_mse ylim_max_mse], ...
     [0.8 0.92 1], 'EdgeColor', 'none','FaceAlpha', 0.5); 
hold on;  

h1_mse = plot(myEpsilon, mse_smooth(1,:), '-.o', 'Color', Color(1,:));
h2_mse = plot(myEpsilon, mse_smooth(2,:), '-.s', 'Color', Color(2,:));
h3_mse = plot(myEpsilon, mse_smooth(3,:), '-.^', 'Color', Color(3,:));
h4_mse = plot(myEpsilon, mse_min_smooth(1,:), '-.p', 'Color', Color(5,:));

xlim([xlim_min, xlim_max]);
ylim([ylim_min_mse, ylim_max_mse]);
set(gca, 'YScale', 'log');

if ylim_max_mse/ylim_min_mse > 100
    yticks([0.1 1 10 100 1000 10000 100000]); 
else
    log_min = floor(log10(ylim_min_mse));
    log_max = ceil(log10(ylim_max_mse));
    ytick_values = 10.^(log_min:log_max);
    yticks(ytick_values);
end 

xlabel('$\epsilon$', 'FontSize', 35);
ylabel('Mean MSE', 'FontSize', 20);
legend([h1_mse h2_mse h3_mse h4_mse], {'RAPPOR','OUE','GRR','Mixture'});
grid on;

text('Parent',ah2(1),'string',{sprintf('Mean over %d features', n_features)},...
    'Position',[xlim_min + 0.5*(xlim_max-xlim_min), ylim_max_mse*0.8, 0]);

annotation('textbox',...
    [0.386 0.839909089283511 0.10950000231266 0.0823863652619449],...
    'String',{'High','Privacy'}, 'EdgeColor','none');
annotation('textbox',...
    [0.549 0.217045452919874 0.106500002241135 0.0823863652619449],...
    'String',{'Low','Privacy'}, 'EdgeColor','none');

filename_mse = sprintf('MSE_mean_%d_features.png', n_features);
print(gcf, filename_mse, '-dpng', '-r600');
fprintf('Saved MSE plot as: %s\n', filename_mse);

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

function smoothed_data = fix_anomalies(data, threshold_factor)
    if nargin < 2
        threshold_factor = 3;
    end
    smoothed_data = data;
    n = length(data);
    if n < 3
        return;
    end
    diffs = abs(diff(data));
    median_diff = median(diffs);
    threshold = threshold_factor * median_diff;
    for i = 2:n-1
        left_jump = abs(data(i) - data(i-1));
        right_jump = abs(data(i+1) - data(i));
        if left_jump > threshold && right_jump > threshold
            expected_value = (data(i-1) + data(i+1)) / 2;
            if abs(data(i) - expected_value) > threshold
                fprintf('Detected anomaly at point %d: %.2f -> %.2f\n', i, data(i), expected_value);
                smoothed_data(i) = expected_value;
            end
        end
    end
end