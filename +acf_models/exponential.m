function acf_filter = exponential(X_grid, Y_grid, corr_lengths)
%EXPONENTIAL Creates a 2D exponential filter for the ACF method.
%   Handles both isotropic (scalar corr_lengths) and anisotropic
%   (2-element vector corr_lengths) cases.

if isscalar(corr_lengths)
    clx = corr_lengths;
    cly = corr_lengths;
else
    clx = corr_lengths(1);
    cly = corr_lengths(2);
end

% The filter corresponds to an Exponential Autocovariance Function
F = exp(-(abs(X_grid)/(clx/2) + abs(Y_grid)/(cly/2)));

% Normalize the filter
acf_filter = F / sum(F(:));
end