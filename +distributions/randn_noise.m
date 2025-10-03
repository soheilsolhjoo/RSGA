function noise = randn_noise(num_points)
%GAUSSIAN Generates uncorrelated noise with a standard normal distribution.
%   noise = GAUSSIAN(num_points) returns an [Ny x Nx] matrix of random
%   numbers drawn from a standard normal (Gaussian) distribution.

noise = randn(num_points(2), num_points(1));

end