function analysis_results = analyze_surface(surface_data, config)
%ANALYZE_SURFACE Calculates statistical and spectral properties of a surface.
%   analysis_results = ANALYZE_SURFACE(surface_data, config) takes the
%   surface data and the configuration struct, performs the analyses
%   specified in config.analysis, and returns the results in a new struct.

disp('--- Starting Surface Analysis ---');

surface_h = surface_data.surface_h;
x_coords = surface_data.x_coords;
y_coords = surface_data.y_coords;
[Ny, Nx] = size(surface_h);

analysis_results = struct();

dx = x_coords(2) - x_coords(1);
dy = y_coords(2) - y_coords(1);
Lx = (Nx-1)*dx; 
Ly = (Ny-1)*dy;
analysis_results.grid.dx = dx;
analysis_results.grid.dy = dy;
analysis_results.grid.L = [Lx, Ly];
analysis_results.grid.Nx = Nx;
analysis_results.grid.Ny = Ny;

if config.analysis.calculate_statistics
    disp('  - Calculating height distribution and statistics...');
    Z = surface_h(:);
    analysis_results.statistics.skewness = mean((Z - mean(Z)).^3) / (std(Z)^3);
    analysis_results.statistics.kurtosis = mean((Z - mean(Z)).^4) / (std(Z)^4);
    
    num_bins = 100;
    [counts, bin_centers] = hist(Z, num_bins);
    analysis_results.hist.counts = counts;
    analysis_results.hist.bin_centers = bin_centers;
end

if config.analysis.calculate_moments
    disp('  - Calculating spectral moments...');
    m0 = mean(surface_h(:).^2);
    analysis_results.moments.m0 = m0;
    analysis_results.rms_height = sqrt(m0);
    
    [dfdx, dfdy] = gradient(surface_h, dx, dy);
    df2 = dfdx.^2 + dfdy.^2;
    m2 = mean(df2(:))/2;
    analysis_results.moments.m2 = m2;
    analysis_results.rms_gradient = sqrt(2*m2);
    analysis_results.rms_gradient_x = sqrt(mean(dfdx(:).^2));
    analysis_results.rms_gradient_y = sqrt(mean(dfdy(:).^2));
    
    [d2fdx2, ~] = gradient(dfdx, dx, dy);
    [~, d2fdy2] = gradient(dfdy, dx, dy);
    d2f2 = d2fdx2.^2 + d2fdy2.^2;
    m4 = mean(d2f2(:));
    analysis_results.moments.m4 = m4;
end

if config.analysis.calculate_acf
    disp('  - Calculating autocovariance functions...');
    acfx = zeros(1, Nx);
    acfy = zeros(1, Ny);
    for i = 1:Ny; cx = xcov(surface_h(i,:), 'coeff'); acfx = acfx + cx(Nx:end); end
    for i = 1:Nx; cy = xcov(surface_h(:,i), 'coeff'); acfy = acfy + cy(Ny:end)'; end
    acfx = acfx / Ny; acfy = acfy / Nx;
    analysis_results.acf.acfx = acfx; analysis_results.acf.acfy = acfy;
    analysis_results.acf.lags_x = linspace(0, Lx, Nx); analysis_results.acf.lags_y = linspace(0, Ly, Ny);
    idx_x = find(acfx < 1/exp(1), 1, 'first'); if isempty(idx_x); idx_x = Nx; end
    analysis_results.acf.corr_length_x = analysis_results.acf.lags_x(idx_x);
    idx_y = find(acfy < 1/exp(1), 1, 'first'); if isempty(idx_y); idx_y = Ny; end
    analysis_results.acf.corr_length_y = analysis_results.acf.lags_y(idx_y);
end

if config.analysis.calculate_psd
    disp('  - Calculating power spectral density...');
    
    fft_surface = fftshift(fft2(surface_h));
    
    % --- PSD Calculation from FFT (Analyzer) ---
    scaling_factor = (Lx * Ly) / ((Nx * Ny)^2 * (4*pi^2));
    psd_2d = scaling_factor * abs(fft_surface).^2;
    analysis_results.psd.psd_2d = psd_2d;
    
    % --- 1D PSD using radialavg (binned by equal radius) ---
    [psd_vals, r_indices] = radialavg(psd_2d, floor(min(Nx,Ny)/2));
    if ~isempty(psd_vals); psd_vals(1) = []; r_indices(1) = []; end
    point_spacing_eff = sqrt(dx * dy);
    q_max = pi / point_spacing_eff;
    analysis_results.psd.radialavg.q = r_indices * q_max;
    analysis_results.psd.radialavg.psd = psd_vals;
    
    % --- 1D PSD using RadiAve (binned by integer pixel distance) ---
    [psd_vals_ra, q_indices_ra] = radiave(psd_2d);
    q_step = 2*pi / sqrt(Lx*Ly);
    analysis_results.psd.radiave.q = q_indices_ra' * q_step;
    analysis_results.psd.radiave.psd = psd_vals_ra;
    
    % --- Slices and Grids (for plotting target PSD) ---
    qx_vec = 2*pi * (-Nx/2 : Nx/2-1) / Lx;
    qy_vec = 2*pi * (-Ny/2 : Ny/2-1) / Ly;
    [Qx_grid, Qy_grid] = meshgrid(qx_vec, qy_vec);
    analysis_results.psd.Qx_grid = Qx_grid;
    analysis_results.psd.Qy_grid = Qy_grid;
    center_x_idx = floor(Nx/2) + 1;
    center_y_idx = floor(Ny/2) + 1;
    analysis_results.psd.psd_slice_x = psd_2d(center_y_idx, :);
    analysis_results.psd.psd_slice_y = psd_2d(:, center_x_idx)';
    analysis_results.psd.qx_vec = qx_vec;
    analysis_results.psd.qy_vec = qy_vec;
    
    % disp('  - Calculating power spectral density...');
    % 
    % % Use the correct scaling factor consistent with the generator's constraints
    % fft_surface = fftshift(fft2(surface_h));
    % scaling_factor = (Lx * Ly) / ((Nx * Ny)^2 * (4*pi^2));
    % psd_2d = scaling_factor * abs(fft_surface).^2;
    % analysis_results.psd.psd_2d = psd_2d;
    % 
    % % Create frequency vectors and grids
    % qx_vec = 2*pi * (-Nx/2 : Nx/2-1) / Lx;
    % qy_vec = 2*pi * (-Ny/2 : Ny/2-1) / Ly;
    % 
    % % --- Create and store the 2D grids needed by report_results.m ---
    % [Qx_grid, Qy_grid] = meshgrid(qx_vec, qy_vec);
    % analysis_results.psd.Qx_grid = Qx_grid;
    % analysis_results.psd.Qy_grid = Qy_grid;
    % 
    % % Extract 1D slices for directional analysis
    % center_x_idx = floor(Nx/2) + 1;
    % center_y_idx = floor(Ny/2) + 1;
    % analysis_results.psd.psd_slice_x = psd_2d(center_y_idx, :);
    % analysis_results.psd.psd_slice_y = psd_2d(:, center_x_idx)';
    % analysis_results.psd.qx_vec = qx_vec;
    % analysis_results.psd.qy_vec = qy_vec;
    % 
    % % Calculate radial average for isotropic visualization
    % [psd_radial, r_indices] = radialavg(psd_2d, floor(min(Nx,Ny)/2));
    % if ~isempty(psd_radial)
    %     psd_radial(1) = []; % Remove DC component
    %     r_indices(1) = [];
    % end
    % 
    % point_spacing_eff = sqrt(dx * dy);
    % q_max = pi / point_spacing_eff;
    % analysis_results.psd.q_radial = r_indices * q_max;
    % analysis_results.psd.psd_radial = psd_radial;
    % 
    % % --- Temporary Test: Calculate PSD using RadiAve.m for comparison ---
    % disp('  - Calculating PSD using RadiAve method for comparison...');
    % psd_radial_RadiAve = radiave(psd_2d);
    % 
    % % The RadiAve function does not return a q-vector, so we create one.
    % % It returns Ny/2 points. We will create a linspace for its q-axis.
    % q_lower_est = 2*pi / max(Lx, Ly);
    % q_upper_est = pi / min(dx, dy);
    % analysis_results.psd.q_radial_RadiAve = linspace(q_lower_est, q_upper_est, length(psd_radial_RadiAve));
    % analysis_results.psd.psd_radial_RadiAve = psd_radial_RadiAve;
end

% --- Spectral Moments Calculation from Numerical PSD ---
if isfield(analysis_results, 'psd') && isfield(analysis_results.psd, 'radialavg')
    disp('  - Calculating spectral moments from numerical PSD...');

    % q_vec = analysis_results.psd.radiave.q;
    % psd_vec = analysis_results.psd.radiave.psd;
    q_vec = analysis_results.psd.radialavg.q;
    psd_vec = analysis_results.psd.radialavg.psd;

    % Ensure no NaNs or Infs are in the data to be integrated
    valid_idx = isfinite(q_vec) & isfinite(psd_vec);
    q_vec = q_vec(valid_idx);
    psd_vec = psd_vec(valid_idx);

    % Define integrands
    integrand_m0 = 2 * pi * q_vec .* psd_vec;
    integrand_m2 = pi * q_vec.^3 .* psd_vec;
    integrand_m4 = 3 * pi / 4 * q_vec.^5 .* psd_vec;

    % Perform numerical integration using the trapezoidal method
    m0_spectral = trapz(q_vec, integrand_m0);
    m2_spectral = trapz(q_vec, integrand_m2);
    m4_spectral = trapz(q_vec, integrand_m4);

    % Store the results
    analysis_results.moments.m0_spectral = m0_spectral;
    analysis_results.moments.m2_spectral = m2_spectral;
    analysis_results.moments.m4_spectral = m4_spectral;
end

if config.analysis.calculate_asperity_params
    disp('  - Calculating multi-asperity parameters...');

    % Ensure moments are calculated first
    if ~exist('m0', 'var') || ~exist('m2', 'var') || ~exist('m4', 'var')
        error('Spectral moments must be calculated before asperity parameters. Set config.analysis.calculate_moments = true;');
    end

    alpha_param = (m0 * m4) / m2^2;

    % Store results in a sub-struct
    analysis_results.asperity.alpha_param = alpha_param;
    analysis_results.asperity.summit_density = m4 / (6*pi*sqrt(3)*m2); % Eta from
    analysis_results.asperity.mean_summit_radius = (3/8)*sqrt(pi/m4); % R from
    analysis_results.asperity.summit_height_std = sqrt(m0 * (1 - 0.8968/alpha_param));
    analysis_results.asperity.mean_summit_height = 4 * sqrt(m0 / (pi * alpha_param));
end

% --- Hybrid Roughness Parameter Calculation ---
if isfield(analysis_results, 'rms_height') && isfield(analysis_results, 'rms_gradient') && isfield(analysis_results, 'acf')
    % Use the arithmetic mean for the effective correlation length
    effective_cl = (analysis_results.acf.corr_length_x + analysis_results.acf.corr_length_y) / 2;
    
    if effective_cl > 0
        analysis_results.hybrid_roughness_param = analysis_results.rms_height * analysis_results.rms_gradient / effective_cl;
    else
        analysis_results.hybrid_roughness_param = inf;
    end
end

disp('--- Surface Analysis Complete ---');
end