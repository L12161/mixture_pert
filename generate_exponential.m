function randnum = generate_exponential(num, ub)
% Generates random numbers following an exponential distribution.
%
% Parameters:
%   num: Number of random values to generate
%   ub: Upper bound for the generated values
%
% Returns:
%   randnum: Array of random integers following an exponential distribution

    % Generate exponential with mean ub/5 (to ensure most values fall within [1, ub])
    scale = ub / 5;
    
    % Generate continuous exponential random variables
    continuous = exprnd(scale, num, 1);
    
    % Discretize and clip to [1, ub] range
    randnum = round(continuous);
    randnum = max(1, min(ub, randnum));
end
