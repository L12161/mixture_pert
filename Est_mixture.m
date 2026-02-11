function est_count = est_mixture(W_list, N_loc, N_user, a_list, b_list, p_list, q_list, alpha_param, count_true)
% Frequency estimator for mixture GRR-BUE mechanism (unbiased)
% Formula: ĉᵢ = (cobs^i - n[π*qᵢ + (1-π)*bᵢ]) / [π(pᵢ - qᵢ) + (1-π)(aᵢ - bᵢ)]
% 
% Inputs:
%   W_list - privacy level assignment for each attribute
%   N_loc - number of attributes
%   N_user - total number of users
%   a_list, b_list - BUE mechanism parameters (size: N_loc)
%   p_list, q_list - GRR mechanism parameters (size: N_loc)
%   alpha_param - mixing probability π (scalar)
%   count_true - true frequency counts
%
% Outputs:
%   est_count - estimated frequency counts (unbiased)

est_count = zeros(N_loc,1);

parfor k = 1:N_loc
    % Get true count for attribute k
    n1 = count_true(k);
    n2 = N_user - n1;
    
    % Get privacy parameters for this attribute
    a = a_list(W_list(k));
    b = b_list(W_list(k));
    p = p_list(W_list(k));
    q = q_list(W_list(k));
    
    % Generate responses using mixture mechanism:
    % With probability alpha_param, user follows GRR
    % With probability (1-alpha_param), user follows BUE
    
    % For users with true value 1:
    rand_vals_1 = rand(n1,1);
    uses_grr_1 = rand_vals_1 < alpha_param;
    
    y1 = zeros(n1,1);
    y1(uses_grr_1) = randsrc(sum(uses_grr_1),1,[[0 1]; [1-p p]]);  % GRR: report 1 with prob p
    y1(~uses_grr_1) = randsrc(sum(~uses_grr_1),1,[[0 1]; [1-a a]]);  % BUE: report 1 with prob a
    
    % For users with true value 0:
    rand_vals_2 = rand(n2,1);
    uses_grr_2 = rand_vals_2 < alpha_param;
    
    y2 = zeros(n2,1);
    y2(uses_grr_2) = randsrc(sum(uses_grr_2),1,[[0 1]; [1-q q]]);  % GRR: report 1 with prob q
    y2(~uses_grr_2) = randsrc(sum(~uses_grr_2),1,[[0 1]; [1-b b]]);  % BUE: report 1 with prob b
    
    % Observed count
    raw_count = sum([y1; y2]);
    
    % Unbiased estimator: ĉᵢ = (cobs^i - n[π*qᵢ + (1-π)*bᵢ]) / [π(pᵢ - qᵢ) + (1-π)(aᵢ - bᵢ)]
    numerator = raw_count - N_user * (alpha_param * q + (1 - alpha_param) * b);
    denominator = alpha_param * (p - q) + (1 - alpha_param) * (a - b);
    
    est_count(k) = numerator / denominator;
end

end