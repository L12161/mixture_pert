function randnum = generate_zipf(num, ub)
% Generates random numbers following a Zipf distribution.
%
% Parameters:
%   num: Number of random values to generate
%   ub: Upper bound for the generated values
%
% Returns:
%   randnum: Array of random integers following a Zipf distribution

    % Zipf with parameter alpha (smaller alpha means more skewed)
    alpha = 1.5;
    
    % Generate numbers from 1 to ub
    x = (1:ub)';
    
    % Calculate Zipf PMF
    pmf = 1 ./ (x .^ alpha);
    
    % Normalize to get probabilities
    pmf = pmf / sum(pmf);
    
    % Generate random indices based on the PMF
    % Use custom sampling function for MATLAB
    randnum = sample_discrete(pmf, num);
end