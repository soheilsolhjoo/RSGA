function [avg_values, q_indices] = radiave(matrix_2d)
% Calculates radial average of a 2D matrix using Euclidean distance.
% This version uses accumarray for efficient averaging.

N = size(matrix_2d, 1);
center = ceil(N / 2);

% Create a grid of pixel indices
[cols, rows] = meshgrid(1:N, 1:N);

% Calculate Euclidean distance and round to nearest integer to define bins
pixel_dist_bins = round(sqrt((rows - center).^2 + (cols - center).^2)) + 1;

% Use accumarray to efficiently calculate the mean of all values in each bin
% It groups all values from matrix_2d based on their bin number in pixel_dist_bins
avg_values = accumarray(pixel_dist_bins(:), matrix_2d(:), [], @mean);

% Trim the output to N/2, as in the original function
num_output_points = ceil(N / 2);
if length(avg_values) > num_output_points
    avg_values = avg_values(1:num_output_points);
else
    % Pad with zeros if there are fewer unique distances than N/2
    avg_values(end+1:num_output_points) = 0;
end

% The q_indices are the bin numbers (0, 1, 2, ...)
q_indices = (0:length(avg_values)-1)';
avg_values = avg_values'; % Return as row vector for consistency
end