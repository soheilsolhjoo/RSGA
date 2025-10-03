function [C0, params] = solve_psd_constraint(config, grid_params)
%SOLVE_PSD_CONSTRAINT Solves for C0 and qr using analytical formulas.
%   [C0, params] = SOLVE_PSD_CONSTRAINT(config, q_grids) calculates the
%   necessary amplitude C0 and/or roll-off qr to satisfy the specified
%   constraint mode.

params = config.psd;
psd_slope = params.psd_slope;
q_lower = grid_params.q_lower;
q_upper = grid_params.q_upper;

params.q_lower = q_lower;
params.q_upper = q_upper;

assert(isscalar(params.corr_lengths), 'For this PSD model, config.psd.corr_lengths must be a scalar.');
cl_eff = params.corr_lengths;
q_rolloff = 2*pi / cl_eff;
params.q_rolloff = q_rolloff;

A0 = 2*pi; A2 = pi;
Ml = q_lower / q_rolloff; Ms = q_upper / q_rolloff;
I0 = ((1 - Ml^2) / 2) - ((1 - Ms^(2 + psd_slope)) / (2 + psd_slope + eps));
I2 = ((1 - Ml^4) / 4) - ((1 - Ms^(4 + psd_slope)) / (4 + psd_slope + eps));

switch params.constraint_mode
    case 'm0'
        sigma_rms = params.target_rms_height;
        C0 = sigma_rms^2 / (A0 * I0 * q_rolloff^2);
    case 'm2'
        g_rms = params.target_rms_gradient;
        C0 = g_rms^2 / 2 / (A2 * I2 * q_rolloff^4);
    case 'm0_and_m2'
        sigma_rms = params.target_rms_height;
        g_rms = params.target_rms_gradient;
        qr_solved = (g_rms / sigma_rms) * sqrt(I0 / I2);
        params.q_rolloff = qr_solved;
        Ml = q_lower / qr_solved;
        Ms = q_upper / qr_solved;
        I0_solved = ((1 - Ml^2) / 2) - ((1 - Ms^(2 + psd_slope)) / (2 + psd_slope + eps));
        C0 = sigma_rms^2 / (A0 * I0_solved * qr_solved^2);
end
end