function psd_2d = simple_rolloff(Qx_grid, Qy_grid, C0, params)
%SIMPLE_ROLLOFF Calculates a 2D PSD for a fractal surface with a roll-off.
%   psd_2d = SIMPLE_ROLLOFF(Qx_grid, Qy_grid, C0, params) returns a 2D PSD
%   based on the provided frequency grids and parameters. The PSD is flat
%   below a roll-off frequency and follows a power-law decay above it.
%   
%
%   Inputs:
%       Qx_grid, Qy_grid - 2D frequency coordinate grids.
%       C0               - The amplitude of the PSD at the roll-off frequency.
%       params           - Struct with fields: .psd_slope, .rolloff_freq.

psd_slope = params.psd_slope;
q_lower = params.q_lower;
q_upper = params.q_upper;
q_rolloff = params.q_rolloff;

q_mag = sqrt(Qx_grid.^2 + Qy_grid.^2);
% --- Replace the q=0 value at the center ---
% Find the indices of the central (DC) component
center_y = floor(size(q_mag, 1) / 2) + 1;
center_x = floor(size(q_mag, 2) / 2) + 1;

% Replace the q=0 value with the value of its immediate neighbor
% This regularizes the q=0 point without affecting the rest of the spectrum.
q_mag(center_y, center_x) = q_mag(center_y, center_x + 1);

psd_2d = zeros(size(Qx_grid));

idx_region2 = (q_mag >= q_lower) & (q_mag <= q_rolloff);
psd_2d(idx_region2) = C0;

idx_region3 = (q_mag > q_rolloff) & (q_mag <= q_upper);
psd_2d(idx_region3) = C0 * (q_mag(idx_region3) / q_rolloff).^psd_slope;

end