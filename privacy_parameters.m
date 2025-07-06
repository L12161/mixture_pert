function [a_values, b_values, alpha_values, p_values, q_values] = privacy_parameters(W, Prob, d)
% ANALYZE_PRIVACY_PARAMETERS Calculate privacy mechanism parameters
%   [a, b, alpha, p, q] = analyze_privacy_parameters(W, Prob, N_user)
%
% Inputs:
%   W - Array of privacy weight levels (e.g., [1 1.2 2])
%   Prob - Probability distribution for each level (e.g., [0.05 0.05 0.9])
%   N_user - Number of users in the dataset
%
% Outputs:
%   a_values - Optimal a parameters for each epsilon level
%   b_values - Optimal b parameters for each epsilon level
%   alpha_values - Alpha parameters for mixture mechanism
%   p_values - p parameters for GRR component
%   q_values - q parameters for GRR component

% Initialize parameters
N_lev = length(W);
N_loc = 100;  % Default location parameter
alpha = N_loc * Prob;

% Load or generate data (placeholder - modify as needed)
% data = generate_exponential(N_user, N_loc);

% Define optimization functions
fun0 = @(x) alpha * ((x(N_lev+1:end) - x(N_lev+1:end).^2) ./ ((x(1:N_lev)-x(N_lev+1:end)).^2)) ...
    + max((1-x(1:N_lev)-x(N_lev+1:end)) ./ (x(1:N_lev)-x(N_lev+1:end)));

fun1 = @(x) alpha * (exp(x) ./ ((exp(x)-1).^2)); % symmetric
fun2 = @(b) alpha * ((b-b.^2) ./ ((0.5-b).^2)) + 1; % a = 0.5

% Epsilon values to evaluate
myEpsilon = 1;
N_epsilon = length(myEpsilon);

% Initialize output variables
a_values = [];
b_values = [];
alpha_values = [];
p_values = [];
q_values = [];

% Main optimization loop
for i = 1:N_epsilon
    epsilon = myEpsilon(i);
    temp = W * epsilon;
    
    % Optimization for mixture mechanism
    % fun3 = @(x)  x(4*N_lev+1)*sum((x(2*N_lev+1 : 3*N_lev).*(1-x(2*N_lev+1 : 3*N_lev)))./(x(2*N_lev+1 : 3*N_lev)-x(3*N_lev+1 : 4*N_lev)).^2) ...
    %     +...
    %     alpha*((( ((1 - x(4*N_lev+1:end)) ).*x(N_lev+1:2*N_lev) - ((1 - x(4*N_lev+1:end))).*x(N_lev+1:2*N_lev).^2 )) ...
    %     ./ ...
    %     ((x(1:N_lev)-x(N_lev+1:2*N_lev)).^2))...
    %     + ...
    %     max( ((1 - x(4*N_lev+1:end)) - (1 - x(4*N_lev+1:end)).*x(1:N_lev) - (1 - x(4*N_lev+1:end)).*x(N_lev+1:2*N_lev)) ...
    %     ./ ...
    %     (x(1:N_lev)-x(N_lev+1:2*N_lev)) );
    fun3 = @(x)  x(4*N_lev+1)*sum( (x(3*N_lev+1 : 4*N_lev).*(1-x(3*N_lev+1 : 4*N_lev)))... 
        ./...
        (x(4*N_lev+1)*(x(2*N_lev+1 : 3*N_lev)-x(3*N_lev+1 : 4*N_lev)).^2 + (1-x(4*N_lev+1))*(x(1:N_lev)-x(N_lev+1:2*N_lev)).^2) ...
        )...
        + ... % sum term for IDUE 
        alpha*((( ((1 - x(4*N_lev+1:end)) ).*x(N_lev+1:2*N_lev) - ((1 - x(4*N_lev+1:end))).*x(N_lev+1:2*N_lev).^2 )) ...
        ./ ...
        (x(4*N_lev+1)*(x(2*N_lev+1 : 3*N_lev)-x(3*N_lev+1 : 4*N_lev)).^2 + (1-x(4*N_lev+1))*(x(1:N_lev)-x(N_lev+1:2*N_lev)).^2))...
        + ... % max term from IDUE
        max( (x(1:N_lev)-x(N_lev+1:2*N_lev)).*((1 - x(4*N_lev+1:end)) - (1 - x(4*N_lev+1:end)).*x(1:N_lev) - (1 - x(4*N_lev+1:end)).*x(N_lev+1:2*N_lev)) ...
        ./ ...
        (x(4*N_lev+1)*(x(2*N_lev+1 : 3*N_lev)-x(3*N_lev+1 : 4*N_lev)).^2 + (1-x(4*N_lev+1))*(x(1:N_lev)-x(N_lev+1:2*N_lev)).^2)) ...
        + ...% max term from GRR
        max( ( x(2*N_lev+1 : 3*N_lev) - x(3*N_lev+1 : 4*N_lev) ).*( x(4*N_lev+1:end) - x(4*N_lev+1:end).*x(2*N_lev+1 : 3*N_lev) -  x(4*N_lev+1:end).*x(3*N_lev+1 : 4*N_lev) )  ...
        ./ ...
        (x(4*N_lev+1)*(x(2*N_lev+1 : 3*N_lev)-x(3*N_lev+1 : 4*N_lev)).^2 + (1-x(4*N_lev+1))*(x(1:N_lev)-x(N_lev+1:2*N_lev)).^2) );
    
    [X, ~] = min_opt3(temp, fun3, d);
    
    % Store results
    a_values = [a_values; X(1:N_lev)];
    b_values = [b_values; X(N_lev+1:2*N_lev)];
    p_values = [p_values; X(2*N_lev+1:3*N_lev)];
    q_values = [q_values; X(3*N_lev+1:4*N_lev)];
    alpha_values = [alpha_values; X(4*N_lev+1:end)];
end
end