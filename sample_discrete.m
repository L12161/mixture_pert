function samples = sample_discrete(pmf, n)
% Samples from a discrete distribution given by PMF
%
% Parameters:
%   pmf: Probability mass function (must sum to 1)
%   n: Number of samples to draw
%
% Returns:
%   samples: n samples from the distribution

    % Create cumulative distribution function
    cdf = cumsum(pmf);
    
    % Generate n uniform random numbers
    u = rand(n, 1);
    
    % Initialize samples
    samples = zeros(n, 1);
    
    % Perform inverse transform sampling
    for i = 1:n
        samples(i) = find(u(i) <= cdf, 1, 'first');
    end
end