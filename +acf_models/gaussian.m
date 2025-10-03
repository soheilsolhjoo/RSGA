function acf_filter = gaussian(X_grid, Y_grid, corr_lengths)
%GAUSSIAN Creates a 2D Gaussian filter for the ACF method.
%   acf_filter = GAUSSIAN(X_grid, Y_grid, corr_lengths) returns a 2D
%   Gaussian filter based on the provided coordinate grids and correlation
%   lengths.
%
%   Handles both isotropic (scalar corr_lengths) and anisotropic
%   (2-element vector corr_lengths) cases.

if isscalar(corr_lengths)
    % Isotropic case
    clx = corr_lengths;
    cly = corr_lengths;
else
    % Anisotropic case
    clx = corr_lengths(1);
    cly = corr_lengths(2);
end

% The filter corresponds to a Gaussian Autocovariance Function.
F = exp(-((X_grid.^2 / (clx^2 / 2)) + (Y_grid.^2 / (cly^2 / 2))));

% Normalize the filter
acf_filter = F / sum(F(:));

end