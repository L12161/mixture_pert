function randnum = generate_uniform(num, ub)
% Generates random numbers following a uniform distribution.
%
% Parameters:
%   num: Number of random values to generate
%   ub: Upper bound for the generated values
%
% Returns:
%   randnum: Array of random integers following a uniform distribution

    % Generate uniform random integers between 1 and ub
    randnum = randi(ub, num, 1);
end