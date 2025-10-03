function y = gaussian(x, params)
%GAUSSIAN Returns the PDF of a Gaussian distribution.
%   y = GAUSSIAN(x, params) calculates the probability density of a
%   Gaussian distribution at the points in x.
%
%   params is a struct containing:
%       .rms_height  (Standard deviation sigma)

mu = 0; % Mean is assumed to be zero
sigma = params.rms_height;
y = normpdf(x, mu, sigma);
end