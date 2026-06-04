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

%%

N_lev = 3;  % the different stages of epsilon values
W = [1 1.2 2]; % these terms are? 
N_loc = 100;    
% used as upper bound when generating synthetic data 
N_user = 100000; 
% no of synthetic users 
Prob = [0.05 0.05 0.9];
alpha = N_loc*Prob; 
% alpha refers to the privacy budget distribution. There are 3 levels of privacy budget here. 5% of them are epsilon, 5% are 1.2 epsilon, 90% are 2epsilon 

W_list = [];
for j = 1:N_lev
    W_list = [W_list; j*ones(alpha(j),1)];
end     
%  W_list contains 5 1's 5 2's and 90 3's after this loop 
W_list = W_list(randperm(N_loc));   
% randomly taken 100 samples from the previous W  
% W_list is used in the actual L1 measurement 
W_tab = tabulate(W_list);
% found no use of W_tab 

W_index = cell(N_lev,1);
for j = 1:N_lev
    W_index{j} = find(W_list==j);
end
% found no use of W_index 

data = generate_powlaw(N_user,N_loc);
% This functions re     ads the unique epsilon values and the proportions of
% these epsilon values for every attribute of a given dataset. 
%domain_size = numel(unique(data));
domain_size = 22980

% ===============  DATASET EVALUATION =======================
% Assume 'data' is your 100000x1 array
unique_vals = unique(data);              % Find all unique values
counts = histc(data, unique_vals);       % Count how many times each value appears
% Mean and standard deviation of the counts
mean_count = mean(counts);
std_count = std(counts);
% Display
fprintf ('Total number of unique counts: %.4f\n', length(unique_vals));
fprintf('Mean of counts: %.4f\n', mean_count);
fprintf('Standard deviation of counts: %.4f\n', std_count);
% ===============  DATASET EVALUATION =======================

% used for actual L1 measurements 

fun0 = @(x) alpha*( (x(N_lev+1:end) - x(N_lev+1:end).^2) ./ ( (x(1:N_lev)-x(N_lev+1:end)).^2 ) )...
    + max( (1-x(1:N_lev)-x(N_lev+1:end)) ./ (x(1:N_lev)-x(N_lev+1:end)) );
% ai -> x(1:N_lev) 
% bi -> x(N_lev+1:end)
% mi -> alpha 

fun1 = @(x) alpha*(exp(x)./((exp(x)-1).^2)); % symmetric
fun2 = @(b) alpha*((b-b.^2)./((0.5-b).^2)) + 1; % a = 0.5

myEpsilon_full = [0.5:0.5:6];
excludeVals = [5];
isExcluded = ismember(myEpsilon_full, excludeVals);
myEpsilon = myEpsilon_full(~isExcluded);
%myEpsilon = [0.5];
% epsilon values starting from 0.5 to 4. Hence, for each of the epsilon
% values, we will generate 5% of epsilon, 5% of 1.2* epsilon and 90% of
% 2*epsilon cases from the user base 
N_epsilon = length(myEpsilon_full);
% epsilon values from 0 to 8 with 0.25 increment 

alpha_values = [];
a_values = [];
b_values = [];
p_values = [];
q_values = [];
for i = 1:N_epsilon
    epsilon = myEpsilon_full(i);
    temp = W*epsilon;
    % here temp will contain 1, 1.2, 2 times epsilon 
    
    % Convert MSE to L1 distance using sqrt approximation
    mse_val = N_loc*exp(epsilon/2)./((exp(epsilon/2)-1).^2);
    result(1,i) = sqrt(mse_val * N_loc);
    a = ones(N_lev,1)*exp(epsilon/2)/(exp(epsilon/2)+1);
    b = 1-a;
    [L1_dist(1,i), ~] = actual_L1(a,b,data,N_loc,W_list);
    
    % result 2 deals with OUE - Convert MSE to L1 distance
    mse_val = N_loc*4*exp(epsilon)./((exp(epsilon)-1).^2)+1;
    result(2,i) = sqrt(mse_val * N_loc);
    a = ones(N_lev,1)*0.5; 
    b = ones(N_lev,1)/(exp(epsilon)+1);
    [L1_dist(2,i), ~] = actual_L1(a,b,data,N_loc,W_list);
    % actual L1 is generating the theoretical L1 values whereas the
    % minopt is generating the empirical values.
    
    % result 3 deals with GRR - Convert MSE to L1 distance
    %mse_val = N_loc * (numel(unique(data)) - 2 + exp(epsilon)) / (exp(epsilon) - 1)^2;
    mse_val = N_loc * (domain_size - 2 + exp(epsilon)) / (exp(epsilon) - 1)^2;
    result(3,i) = sqrt(mse_val * N_loc);

    % Convert optimization results from MSE to L1 distance
    [X, mse_result] = min_opt0(temp,fun0);     % get MSE result
    result_min(1,i) = sqrt(mse_result * N_loc);  % convert to L1 distance
    Xmin0(:,i) = X;  
    a = X(1:N_lev);
    b = X((N_lev+1):end);
    % 3 values from a, 3 from b 
    [L1_dist_min(1,i), ~] = actual_L1(a,b,data,N_loc,W_list);  

    % fun3 corrected with claude AI. The linear sum of variance was kinda correct. 
    fun3 = @(x) x(4*N_lev+1)*sum((x(2*N_lev+1:3*N_lev).*(1-x(2*N_lev+1:3*N_lev)))./(x(2*N_lev+1:3*N_lev)-x(3*N_lev+1:4*N_lev)).^2) ...
        + sum(alpha(:) .* (1-x(4*N_lev+1)) .* x(N_lev+1:2*N_lev).*(1-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)).^2) ...
        + max((1-x(4*N_lev+1)) .* (1-x(1:N_lev)-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)));
    
    % paper notation -> MATLAB Variable name
    %       ai ->  x(1 : N_lev) 
    % bi ->  x(N_lev+1 : 2*N_lev)
    % p -> x(2*N_lev+1 : 3*N_lev)
    % q -> x(3*N_lev+1 : 4*N_lev)
    % Alpha -> x(4*N_lev+1)
    % mi -> alpha 

    % Convert optimization results from MSE to L1 distance
    [X, mse_result] = min_opt3(temp,fun3, domain_size);     
    result_min(4,i) = sqrt(mse_result * N_loc);  % convert to L1 distance
    
    % Extract parameters for mixture mechanism
    a = X(1:N_lev);
    b = X(N_lev+1:2*N_lev);
    p = X(2*N_lev+1:3*N_lev);
    q = X(3*N_lev+1:4*N_lev);
    alpha_param = X(4*N_lev+1);
    
    % Store values for analysis
    a_values = [a_values;a];
    b_values = [b_values;b];
    alpha_values = [alpha_values;alpha_param];
    p_values = [p_values; p];
    q_values = [q_values; q];
    
    % Compute L1 distance using UNBIASED MIXTURE ESTIMATOR
    [L1_dist_min(4,i), ~] = actual_L1_mixture(a,b,p,q,alpha_param,data,N_loc,W_list);
    
    [X, mse_result] = min_opt1(temp,fun1); 
    result_min(2,i) = sqrt(mse_result * N_loc);  % convert to L1 distance
    Xmin1(:,i) = X;
    a = exp(X)./(1+exp(X)); 
    b = 1./(1+exp(X));
    [L1_dist_min(2,i), ~] = actual_L1(a,b,data,N_loc,W_list);
    
    [X, mse_result] = min_opt2(temp,fun2); 
    result_min(3,i) = sqrt(mse_result * N_loc);  % convert to L1 distance
    Xmin2(:,i) = X;
    a = 0.5*ones(N_lev,1); 
    b = X;
    [L1_dist_min(3,i), ~] = actual_L1(a,b,data,N_loc,W_list);
  
end

% =============== COMPUTE AND PRINT MSE ===================
% MSE is the squared L1 distance
mse_result = result.^2;
mse_L1_dist = L1_dist.^2;
mse_result_min = result_min.^2;
mse_L1_dist_min = L1_dist_min.^2;

fprintf('\n============= THEORETICAL MSE (Variance) =============\n');
fprintf('Epsilon\t\tRAPPOR\t\tOUE\t\tGRR\t\tMixture\n');
for i = 1:N_epsilon
    fprintf('%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n', myEpsilon_full(i), ...
        mse_result(1,i), mse_result(2,i), mse_result(3,i), mse_result_min(4,i));
end



%% ===================== L1 DISTANCE PLOT =====================
figure('Position',[100,100,1000,1000]);
ah1 = TightPlots(1, 1, 800,[10 7],[80,80],[50,50,60,40],[70,10],'pixels');
axes(ah1(1));

xlim_min = 0.25;  xlim_max = 8; 
ylim_min = 10; ylim_max = 800;

% Draw background first (for x < 2)
fill([xlim_min 2 2 xlim_min], [ylim_min ylim_min ylim_max ylim_max], ...
     [0.8 0.92 1], 'EdgeColor', 'none','FaceAlpha', 0.5); 
hold on;  % Keep plot active for the next lines

%% Plot continuous lines through all points (including excluded)
h1=plot(myEpsilon_full, result(1,:), '-.', 'Color', Color(1,:), 'LineWidth', 1.3);
h2=plot(myEpsilon_full, result(2,:), '-.', 'Color', Color(2,:), 'LineWidth', 1.3);
h3=plot(myEpsilon_full, result(3,:), '-.', 'Color', Color(3,:), 'LineWidth', 1.3);
h4=plot(myEpsilon_full, result_min(4,:), '-.', 'Color', Color(5,:), 'LineWidth', 1.3);

%% Overlay markers for included points (filled markers)
plot(myEpsilon, result(1,~isExcluded), 'o', 'Color', Color(1,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(1,:), 'LineStyle', 'none');
plot(myEpsilon, result(2,~isExcluded), 's', 'Color', Color(2,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(2,:), 'LineStyle', 'none');
plot(myEpsilon, result(3,~isExcluded), '^', 'Color', Color(3,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(3,:), 'LineStyle', 'none');
plot(myEpsilon, result_min(4,~isExcluded), 'p', 'Color', Color(5,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(5,:), 'LineStyle', 'none');

%% Overlay markers for excluded points (hollow markers)
plot(myEpsilon_full(isExcluded), result(1,isExcluded), 'o', 'Color', Color(1,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(1,:), 'LineStyle', 'none');
plot(myEpsilon_full(isExcluded), result(2,isExcluded), 's', 'Color', Color(2,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(2,:), 'LineStyle', 'none');
plot(myEpsilon_full(isExcluded), result(3,isExcluded), '^', 'Color', Color(3,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(3,:), 'LineStyle', 'none');
plot(myEpsilon_full(isExcluded), result_min(4,isExcluded), 'p', 'Color', Color(5,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(5,:), 'LineStyle', 'none');

% Axes and plot settings
xlim([xlim_min, xlim_max]);
ylim([ylim_min, ylim_max]);
set(gca, 'YScale', 'log');
yticks([0 10 100 1000]); 
xlabel('$\epsilon$', 'FontSize', 37, 'FontWeight', 'bold');
ylabel('L1 Distance', 'FontSize', 24, 'FontWeight', 'bold');
legend([h1 h2 h3 h4], {'RAPPOR','OUE','GRR','Mixture'}, 'FontSize', 16);
grid on;

% Add annotation
text('Parent',ah1(1),'string',{'Uniform Distribution'},...
    'Position',[0.613897420582564 18.4075798429304 0], 'FontSize', 16);

annotation('textbox',...
    [0.386 0.839909089283511 0.10950000231266 0.0823863652619449],...
    'String',{'High','Privacy'}, 'EdgeColor','none', 'FontSize', 16);
annotation('textbox',...
    [0.549 0.217045452919874 0.106500002241135 0.0823863652619449],...
    'String',{'Low','Privacy'}, 'EdgeColor','none', 'FontSize', 16);

% Save figure as high-quality PNG for paper
print(gcf, 'L1_Distance_Comparison.png', '-dpng', '-r600');

%% ===================== MSE PLOT =====================
figure('Position',[100,100,1000,1000]);
ah2 = TightPlots(1, 1, 800,[10 7],[80,80],[50,50,60,40],[70,10],'pixels');
axes(ah2(1));

xlim_min = 0.25;  xlim_max = 8; 

% Calculate adaptive y-limits based on actual MSE data
all_mse_data = [mse_result(1,:), mse_result(2,:), mse_result(3,:), mse_result_min(4,:)];
mse_min_y = min(all_mse_data);
mse_max_y = max(all_mse_data);

% Add some padding (10% on each side)
mse_y_range = mse_max_y - mse_min_y;
mse_ylim_min = max(1, mse_min_y - 0.1*mse_y_range);
mse_ylim_max = mse_max_y + 0.1*mse_y_range;

% Draw background first (for x < 2)
fill([xlim_min 2 2 xlim_min], [mse_ylim_min mse_ylim_min mse_ylim_max mse_ylim_max], ...
     [0.8 0.92 1], 'EdgeColor', 'none','FaceAlpha', 0.5); 
hold on;

%% Plot continuous lines through all points (including excluded)
hm1=plot(myEpsilon_full, mse_result(1,:), '-.', 'Color', Color(1,:), 'LineWidth', 1.3);
hm2=plot(myEpsilon_full, mse_result(2,:), '-.', 'Color', Color(2,:), 'LineWidth', 1.3);
hm3=plot(myEpsilon_full, mse_result(3,:), '-.', 'Color', Color(3,:), 'LineWidth', 1.3);
hm4=plot(myEpsilon_full, mse_result_min(4,:), '-.', 'Color', Color(5,:), 'LineWidth', 1.3);

%% Overlay markers for included points (filled markers)
plot(myEpsilon, mse_result(1,~isExcluded), 'o', 'Color', Color(1,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(1,:), 'LineStyle', 'none');
plot(myEpsilon, mse_result(2,~isExcluded), 's', 'Color', Color(2,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(2,:), 'LineStyle', 'none');
plot(myEpsilon, mse_result(3,~isExcluded), '^', 'Color', Color(3,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(3,:), 'LineStyle', 'none');
plot(myEpsilon, mse_result_min(4,~isExcluded), 'p', 'Color', Color(5,:), 'MarkerSize', 6, 'MarkerFaceColor', Color(5,:), 'LineStyle', 'none');

%% Overlay markers for excluded points (hollow markers)
plot(myEpsilon_full(isExcluded), mse_result(1,isExcluded), 'o', 'Color', Color(1,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(1,:), 'LineStyle', 'none');
plot(myEpsilon_full(isExcluded), mse_result(2,isExcluded), 's', 'Color', Color(2,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(2,:), 'LineStyle', 'none');
plot(myEpsilon_full(isExcluded), mse_result(3,isExcluded), '^', 'Color', Color(3,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(3,:), 'LineStyle', 'none');
plot(myEpsilon_full(isExcluded), mse_result_min(4,isExcluded), 'p', 'Color', Color(5,:), 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Color(5,:), 'LineStyle', 'none');

% Axes and plot settings
xlim([xlim_min, xlim_max]);
ylim([mse_ylim_min, mse_ylim_max]);
set(gca, 'YScale', 'log');
xlabel('$\epsilon$', 'FontSize', 37, 'FontWeight', 'bold');
ylabel('MSE', 'FontSize', 24, 'FontWeight', 'bold');
legend([hm1 hm2 hm3 hm4], {'RAPPOR','OUE','GRR','Mixture'}, 'FontSize', 16);
grid on;

% Add annotation
text('Parent',ah2(1),'string',{'Uniform Distribution'},...
    'Position',[0.613897420582564 18.4075798429304 0], 'FontSize', 16);

annotation('textbox',...
    [0.386 0.839909089283511 0.10950000231266 0.0823863652619449],...
    'String',{'High','Privacy'}, 'EdgeColor','none', 'FontSize', 16);
annotation('textbox',...
    [0.549 0.217045452919874 0.106500002241135 0.0823863652619449],...
    'String',{'Low','Privacy'}, 'EdgeColor','none', 'FontSize', 16);

% Save figure as high-quality PNG for paper
print(gcf, 'MSE_Comparison.png', '-dpng', '-r600');

% Function to compute L1 distance instead of MSE
function [L1_dist, RE] = actual_L1(a_list, b_list, data, N_loc, W_list)
%w list will contain the priv lvl for the users. There should be 90% lvl3 s
%5% lvl 2s and 5% lvl 1s 
repeat = 10;
N_user = length(data);

tab_true = tabulate(data);
% creates a table where the percentage of every entry showing up is
% presented 
count_true = tab_true(:,2);
% count_true holds the true frequency of occurrence of the different classes
L1_dist = zeros(1,repeat);
RE = zeros(N_loc,repeat);

parfor i = 1:repeat
    count_est = est(W_list,N_loc,N_user,a_list,b_list,count_true);
    L1_dist(i) = sum(abs(count_est-count_true));  % L1 distance between estimated and true counts
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

% =====================================================================
% NEW FUNCTION: Unbiased estimator for mixture GRR-BUE mechanism
% =====================================================================
function [L1_dist, RE] = actual_L1_mixture(a_list, b_list, p_list, q_list, alpha_param, data, N_loc, W_list)
% Unbiased estimator for mixture GRR-BUE mechanism
% Formula: ĉᵢ = (cobs^i - n[π*qᵢ + (1-π)*bᵢ]) / [π(pᵢ - qᵢ) + (1-π)(aᵢ - bᵢ)]
% 
% Inputs:
%   a_list, b_list - BUE parameters (N_loc x 1)
%   p_list, q_list - GRR parameters (N_loc x 1)
%   alpha_param - mixing probability π (scalar)
%   data - user data
%   N_loc - number of attributes
%   W_list - privacy level assignment
%
% Outputs:
%   L1_dist - L1 distance (mean over repeats)
%   RE - relative error

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
% Simulated frequency estimation using the unbiased mixture estimator

est_count = zeros(N_loc, 1);

parfor k = 1:N_loc
    n1 = count_true(k);
    n2 = N_user - n1;
    
    % Get parameters for attribute k
    a = a_list(W_list(k));
    b = b_list(W_list(k));
    p = p_list(W_list(k));
    q = q_list(W_list(k));
    
    % Generate responses from mixture mechanism
    % For users with true value 1:
    rand_vals_1 = rand(n1, 1);
    uses_grr_1 = rand_vals_1 < alpha_param;
    
    y1 = zeros(n1, 1);
    y1(uses_grr_1) = randsrc(sum(uses_grr_1), 1, [[0 1]; [1-p p]]);  % GRR
    y1(~uses_grr_1) = randsrc(sum(~uses_grr_1), 1, [[0 1]; [1-a a]]);  % BUE
    
    % For users with true value 0:
    rand_vals_2 = rand(n2, 1);
    uses_grr_2 = rand_vals_2 < alpha_param;
    
    y2 = zeros(n2, 1);
    y2(uses_grr_2) = randsrc(sum(uses_grr_2), 1, [[0 1]; [1-q q]]);  % GRR
    y2(~uses_grr_2) = randsrc(sum(~uses_grr_2), 1, [[0 1]; [1-b b]]);  % BUE
    
    % Observed count
    raw_count = sum([y1; y2]);
    
    % UNBIASED ESTIMATOR: ĉᵢ = (cobs^i - n[π*qᵢ + (1-π)*bᵢ]) / [π(pᵢ - qᵢ) + (1-π)(aᵢ - bᵢ)]
    numerator = raw_count - N_user * (alpha_param * q + (1 - alpha_param) * b);
    denominator = alpha_param * (p - q) + (1 - alpha_param) * (a - b);
    
    est_count(k) = numerator / denominator;
end

end