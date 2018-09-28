function [normed_matrix] = normalize(matrix_to_normalize, wantWhat)
%NORMALIZE Summary of this function goes here
%   Detailed explanation goes here
mean_forEach = mean(matrix_to_normalize, 2);
centered_matrix = bsxfun(@minus, matrix_to_normalize, mean_forEach);
centered_norm_matrix = centered_matrix./std(centered_matrix, 0, 2);
if strcmp(wantWhat, 'centered')
    normed_matrix = centered_norm_matrix;
else
    normed_matrix = bsxfun(@plus, centered_norm_matrix, mean_forEach);
end
end

