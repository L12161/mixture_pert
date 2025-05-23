function randnum = generate_gaussian(num, ub)
% Generates random numbers following a discretized Gaussian distribution.
%
% Parameters:
%   num: Number of random values to generate
%   ub: Upper bound for the generated values
%
% Returns:
%   randnum: Array of random integers following a Gaussian distribution

    % Generate Gaussian centered at ub/2 with standard deviation ub/6
    % This ensures most values fall within [1, ub]
    mean_val = ub / 2;
    std_val = ub / 6;
    
    % Generate continuous Gaussian random variables
    continuous = normrnd(mean_val, std_val, num, 1);
    
    % Discretize and clip to [1, ub] range
    randnum = round(continuous);
    randnum = max(1, min(ub, randnum));
end