function [unique_entries, proportions] = reading_distribution(dataset)
% ANALYZE_COLUMN_ENTRIES Returns unique values and their proportions for each column
%   [unique_entries, proportions] = analyze_column_entries(dataset)
%
% Input:
%   dataset - Table or matrix containing the data to analyze
%
% Outputs:
%   unique_entries - Cell array where each cell contains unique values for a column
%   proportions   - Cell array where each cell contains corresponding proportions

% Convert to table if input is matrix
if ~istable(dataset)
    dataset = array2table(dataset);
end

num_cols = width(dataset);
unique_entries = cell(1, num_cols);
proportions = cell(1, num_cols);

for col = 1:num_cols
    current_col = dataset{:, col};
    
    % Handle both numeric and categorical data
    if iscategorical(current_col)
        current_col = cellstr(current_col);
    end
    
    % Get unique entries and counts
    [unique_vals, ~, idx] = unique(current_col, 'stable');
    counts = accumarray(idx, 1);
    
    % Calculate proportions
    total = sum(counts);
    props = counts / total;
    
    % Store results
    unique_entries{col} = unique_vals;
    proportions{col} = props;
end

end