function psd_2d = k_correlation(Qx_grid, Qy_grid, ~, params)
%PALASANTZAS Calculates a 2D PSD based on the Palasantzas model.
%   The amplitude is determined by the h0 and cl parameters.

h0 = params.target_rms_height;
cl = params.corr_lengths;
psd_slope = params.psd_slope;
q_mag = sqrt(Qx_grid.^2 + Qy_grid.^2);

% Ensure correlation length is a scalar for this isotropic model
if ~isscalar(cl)
    cl = sqrt(cl(1)*cl(2));
    warning('Palasantzas model is isotropic; using geometric mean of correlation lengths.');
end

% The term inside the power can become negative/complex at high q
term = 1 - (cl^2 / (psd_slope + 2)) * q_mag.^2;

% Initialize psd_2d with zeros
psd_2d = zeros(size(q_mag));

% Only calculate PSD where the term is positive
valid_idx = term > 0;
psd_2d(valid_idx) = ((cl*h0)^2 / (2*pi)) * term(valid_idx).^(psd_slope/2);

end