function report_results(surface_data, analysis_results, config)
%REPORT_RESULTS Displays plots and a summary of the analysis.
%   REPORT_RESULTS(surface_data, analysis_results, config) takes the data
%   from generation and analysis and produces visualizations and a text
%   summary based on the switches in the config struct.

disp('--- Reporting Results ---');

% =========================================================================
% --- Plotting ---
% =========================================================================
if isfield(config.plotting, 'show_surface_3d') && config.plotting.show_surface_3d
    figure('Name', 'Generated Surface');
    surf(surface_data.X_grid, surface_data.Y_grid, surface_data.surface_h, 'EdgeColor', 'none');
    axis equal; daspect([1 1 1]);
    title('Generated Surface');
    xlabel('X (grid units)'); ylabel('Y (grid units)'); zlabel('Height (arbitrary units)');
    colorbar; view(3);
end

if isfield(config.plotting, 'show_height_histogram') && config.plotting.show_height_histogram
    figure('Name', 'Height Distribution');
    hold on;

    bar(analysis_results.hist.bin_centers, analysis_results.hist.counts, 'BarWidth', 1, 'FaceAlpha', 0.6);

    Z = surface_data.surface_h(:);
    x_vals = analysis_results.hist.bin_centers;
    bin_width = x_vals(2) - x_vals(1);
    hist_area = numel(Z) * bin_width;

    y_fit = normpdf(x_vals, mean(Z), std(Z));
    plot(x_vals, y_fit * hist_area, 'r--', 'LineWidth', 2);

    target_pdf_handle = [];
    if strcmp(config.generation.method, 'ACF')
        % full_func_name_str = func2str(config.acf.initial_noise_dist);
        % name_parts = strsplit(full_func_name_str, '.');
        target_pdf_handle = str2func(['dist_pdf.' 'gaussian']);
        target_params.rms_height = config.acf.target_rms_height;
    elseif strcmp(config.generation.method, 'Iterative')
        full_func_name_str = func2str(config.iterative.hpd_model);
        name_parts = strsplit(full_func_name_str, '.');
        target_pdf_handle = str2func(['dist_pdf.' name_parts{end}]);
        target_params = config.iterative.hpd_params;
    end

    if ~isempty(target_pdf_handle)
        y_target = target_pdf_handle(x_vals, target_params);
        plot(x_vals, y_target * hist_area, 'k-', 'LineWidth', 2);
        legend('Generated Data', 'Gaussian Fit', 'Target Model');
    else
        legend('Generated Data', 'Gaussian Fit');
    end

    title('Height Distribution Histogram');
    xlabel('Height (arbitrary units)'); ylabel('Counts');
    grid on;
    hold off;
end

if isfield(config.plotting, 'show_acf_plot') && config.plotting.show_acf_plot
    figure('Name', 'Autocovariance Functions');
    min_y = min([analysis_results.acf.acfx, analysis_results.acf.acfy]);
    y_limits = [min_y - 0.1, 1.1];
    subplot(1, 2, 1);
    plot(analysis_results.acf.lags_x, analysis_results.acf.acfx, 'b-');
    hold on;
    title('ACF in X-Direction'); xlabel('Lag (grid units)'); ylabel('Normalized ACF');
    grid on; ylim(y_limits);
    line([analysis_results.acf.corr_length_x, analysis_results.acf.corr_length_x], ylim, 'Color', 'r', 'LineStyle', '--');
    legend('ACF_x', ['CL_x = ' num2str(analysis_results.acf.corr_length_x, '%.2f')], 'Location', 'northeast');
    hold off;
    subplot(1, 2, 2);
    plot(analysis_results.acf.lags_y, analysis_results.acf.acfy, 'b-');
    hold on;
    title('ACF in Y-Direction'); xlabel('Lag (grid units)'); ylabel('Normalized ACF');
    grid on; ylim(y_limits);
    line([analysis_results.acf.corr_length_y, analysis_results.acf.corr_length_y], ylim, 'Color', 'r', 'LineStyle', '--');
    legend('ACF_y', ['CL_y = ' num2str(analysis_results.acf.corr_length_y, '%.2f')], 'Location', 'northeast');
    hold off;
end

if isfield(config.plotting, 'show_psd_plot') && config.plotting.show_psd_plot
    figure('Name', 'Power Spectral Density');
    % hold on;

    % Plot the ACTUAL PSD from the analyzed surface
    if config.grid.is_anisotropic
        title('Anisotropic Power Spectral Density (Slices)');
        idx_x = analysis_results.psd.psd_slice_x > 0;
        loglog(analysis_results.psd.qx_vec(idx_x), analysis_results.psd.psd_slice_x(idx_x), 'b-o', 'MarkerSize', 4, 'DisplayName', 'Actual (Slice along q_x)');
        hold on;
        idx_y = analysis_results.psd.psd_slice_y > 0;
        loglog(analysis_results.psd.qy_vec(idx_y), analysis_results.psd.psd_slice_y(idx_y), 'r-x', 'MarkerSize', 4, 'DisplayName', 'Actual (Slice along q_y)');
    else
        title('Isotropic Power Spectral Density');
        idx1 = analysis_results.psd.radialavg.psd > 0;
        loglog(analysis_results.psd.radialavg.q(idx1), analysis_results.psd.radialavg.psd(idx1), 'b-o', 'DisplayName', 'Actual (Equal-Width Bins)');
        hold on;
        idx2 = analysis_results.psd.radiave.psd > 0;
        loglog(analysis_results.psd.radiave.q(idx2), analysis_results.psd.radiave.psd(idx2), 'g--x', 'DisplayName', 'Actual (Pixel Distance Bins)');
    end

    % --- Plot the TARGET PSD by calling the selected model function ---
    if ~isempty(config) && (strcmp(config.generation.method, 'PSD_FFT') | strcmp(config.generation.method, 'PSD_SUM')) && isfield(surface_data, 'generation_params')

        % 1. Get the model handle from the config used for generation
        model_handle = config.psd.model;

        % 2. Reconstruct the parameters that the model needs
        params = config.psd;
        params.q_rolloff = surface_data.generation_params.q_rolloff;
        params.q_lower = surface_data.generation_params.q_lower;
        params.q_upper = surface_data.generation_params.q_upper;
        C0 = surface_data.generation_params.C0;

        % 3. Use the q-vector from the analysis for the x-axis
        q_vec_target = analysis_results.psd.radialavg.q;

        % 4. Call the model function to get the 1D theoretical PSD values.
        % We pass the 1D q-vector as a dummy Qx_grid and zeros for Qy_grid.
        % The model function will calculate q_mag = sqrt(q_vec_target.^2 + 0.^2),
        % correctly evaluating the isotropic model along the 1D vector.
        psd_target_1d = model_handle(q_vec_target, zeros(size(q_vec_target)), C0, params);

        % 5. Plot the theoretical target curve
        idx_plot = psd_target_1d > 0;
        loglog(q_vec_target(idx_plot), psd_target_1d(idx_plot), 'k--', 'LineWidth', 2, 'DisplayName', 'Target Model');
    end

    xlabel('Wavevector q (rad/grid unit)');
    ylabel('PSD (arbitrary units)');
    grid on;
    % legend(legend_handles, 'Location', 'southwest');
    legend()
    hold off;
end

% =========================================================================
% --- Command Window Summary ---
% =========================================================================
% --- Open file for writing if requested ---
fileID = 1; % Default: print to command window
if isfield(config.output, 'save_summary_txt') && config.output.save_summary_txt

    output_folder = config.output.folder;
    if ~exist(output_folder, 'dir'); mkdir(output_folder); end

    prefix = config.output.filename_prefix;
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = [prefix, '_', timestamp, '.txt'];
    full_filepath = fullfile(output_folder, filename);

    fileID = fopen(full_filepath, 'w');
    if fileID == -1
        warning('Could not open file %s for writing. Summary will only be printed to the command window.', full_filepath);
        fileID = 1; % Revert to command window if file opening fails
    end
end

fprintf(fileID, '\n================== SURFACE ANALYSIS SUMMARY ==================\n');
fprintf(fileID, 'Grid Properties:\n');
fprintf(fileID, '  - Size: %d x %d points\n', analysis_results.grid.Nx, analysis_results.grid.Ny);
fprintf(fileID, '  - Spacing (dx, dy): (%.2f, %.2f)\n', analysis_results.grid.dx, analysis_results.grid.dy);

if isfield(analysis_results, 'statistics')
    fprintf(fileID, '\nStatistical Properties:\n');
    actual_rms_h = analysis_results.rms_height;

    % --- Check if config exists before showing target ---
    if ~isempty(config)
        is_h_constrained = false;
        target_rms_h = 0;
        if strcmp(config.generation.method, 'ACF'); is_h_constrained = true; target_rms_h = config.acf.target_rms_height;
        elseif strcmp(config.generation.method, 'PSD') && (strcmp(config.psd.constraint_mode, 'm0') || strcmp(config.psd.constraint_mode, 'm0_and_m2')); is_h_constrained = true; target_rms_h = config.psd.target_rms_height;
        end
        if is_h_constrained; err_h = abs(actual_rms_h - target_rms_h) / target_rms_h * 100; fprintf(fileID, '  - RMS Height (Sq): %.4f (Target: %.4f, Error: %.2f%%)\n', actual_rms_h, target_rms_h, err_h);
        else; fprintf(fileID, '  - RMS Height (Sq): %.4f\n', actual_rms_h); end
    else
        fprintf(fileID, '  - RMS Height (Sq): %.4f\n', actual_rms_h);
    end
    fprintf(fileID, '  - Skewness (Ssk): %.4f (Gaussian = 0)\n', analysis_results.statistics.skewness);
    fprintf(fileID, '  - Kurtosis (Sku): %.4f (Gaussian = 3)\n', analysis_results.statistics.kurtosis);

    if isfield(analysis_results, 'hybrid_roughness_param')
        fprintf(fileID, '  - Hybrid Parameter (rho):       %.4f\n', analysis_results.hybrid_roughness_param);
    end
end

if isfield(analysis_results, 'moments')
    fprintf(fileID, '\nGradient Properties:\n');
    actual_rms_g = analysis_results.rms_gradient;

    % --- Check if config exists before showing target ---
    if ~isempty(config)
        is_g_constrained = strcmp(config.generation.method, 'PSD') && (strcmp(config.psd.constraint_mode, 'm2') || strcmp(config.psd.constraint_mode, 'm0_and_m2'));
        if is_g_constrained; target_rms_g = config.psd.target_rms_gradient; err_g = abs(actual_rms_g - target_rms_g) / target_rms_g * 100; fprintf(fileID, '  - Eff. RMS Gradient (Sdq): %.4f (Target: %.4f, Error: %.2f%%)\n', actual_rms_g, target_rms_g, err_g);
        else; fprintf(fileID, '  - Eff. RMS Gradient (Sdq): %.4f\n', actual_rms_g); end
    else
        fprintf(fileID, '  - Eff. RMS Gradient (Sdq): %.4f\n', actual_rms_g);
    end
    fprintf(fileID, '  - X-dir RMS Gradient:       %.4f\n', analysis_results.rms_gradient_x);
    fprintf(fileID, '  - Y-dir RMS Gradient:       %.4f\n', analysis_results.rms_gradient_y);
end

if isfield(analysis_results, 'acf')
    fprintf(fileID, '\nCorrelation Properties:\n');
    actual_clx = analysis_results.acf.corr_length_x;
    actual_cly = analysis_results.acf.corr_length_y;

    % --- Check if config exists before showing target ---
    if ~isempty(config)
        target_cl = [];
        if strcmp(config.generation.method, 'ACF'); target_cl = config.acf.corr_lengths;
        elseif strcmp(config.generation.method, 'PSD'); target_cl = config.psd.corr_lengths; end
        if ~isempty(target_cl)
            if isscalar(target_cl); target_clx = target_cl; target_cly = target_cl;
            else; target_clx = target_cl(1); target_cly = target_cl(2); end
            err_clx = abs(actual_clx - target_clx) / target_clx * 100;
            err_cly = abs(actual_cly - target_cly) / target_cly * 100;
            fprintf(fileID, '  - Correlation Length X (Clx): %.2f (Target: %.2f, Error: %.2f%%)\n', actual_clx, target_clx, err_clx);
            fprintf(fileID, '  - Correlation Length Y (Cly): %.2f (Target: %.2f, Error: %.2f%%)\n', actual_cly, target_cly, err_cly);
        else
            fprintf(fileID, '  - Correlation Length X (Clx): %.2f\n', actual_clx);
            fprintf(fileID, '  - Correlation Length Y (Cly): %.2f\n', actual_cly);
        end
    else
        fprintf(fileID, '  - Correlation Length X (Clx): %.2f\n', actual_clx);
        fprintf(fileID, '  - Correlation Length Y (Cly): %.2f\n', actual_cly);
    end
end

if ~isempty(config) && strcmp(config.generation.method, 'PSD') && isfield(surface_data, 'generation_params')
    fprintf(fileID, '\nCharacteristic Wavenumbers (from Generation):\n');
    fprintf(fileID, '  - C0 (Solved Amplitude):  %.4e\n', surface_data.generation_params.C0);
    fprintf(fileID, '  - Lower Cutoff (q_L):   %.4f\n', surface_data.generation_params.q_lower);
    fprintf(fileID, '  - Roll-off (q_r):       %.4f\n', surface_data.generation_params.q_rolloff);
    fprintf(fileID, '  - Upper Cutoff (q_H):   %.4f\n', surface_data.generation_params.q_upper);
end

% if isfield(analysis_results, 'moments')
%     fprintf(fileID, '\nSpectral Moments:\n');
%     fprintf(fileID, '  - m0: %.4e\n', analysis_results.moments.m0);
%     fprintf(fileID, '  - m2: %.4e\n', analysis_results.moments.m2);
%     fprintf(fileID, '  - m4: %.4e\n', analysis_results.moments.m4);
% end

if isfield(analysis_results, 'moments')
    fprintf(fileID, '\nSpectral Moments:\n');
    m = analysis_results.moments;

    err = abs(m.m0 - m.m0_spectral) / m.m0_spectral * 100;
    fprintf(fileID, '  - m0 (Actual):   %.4e (from PSD: %.4e, Diff: %.2f%%)\n', m.m0, m.m0_spectral, err);
    err = abs(m.m2 - m.m2_spectral) / m.m2_spectral * 100;
    fprintf(fileID, '  - m2 (Actual):   %.4e (from PSD: %.4e, Diff: %.2f%%)\n', m.m2, m.m2_spectral, err);
    err = abs(m.m4 - m.m4_spectral) / m.m4_spectral * 100;
    fprintf(fileID, '  - m4 (Actual):   %.4e (from PSD: %.4e, Diff: %.2f%%)\n', m.m4, m.m4_spectral, err);
end

if isfield(analysis_results, 'asperity')
    fprintf(fileID, '\nMulti-Asperity Parameters (GW Model):\n');
    fprintf(fileID, '  - Bandwidth Parameter (alpha):  %.4f\n', analysis_results.asperity.alpha_param);
    fprintf(fileID, '  - Summit Density (eta):         %.4e\n', analysis_results.asperity.summit_density);
    fprintf(fileID, '  - Mean Summit Radius (R):       %.4f\n', analysis_results.asperity.mean_summit_radius);
    fprintf(fileID, '  - Summit Height Std (sigma_s):  %.4f\n', analysis_results.asperity.summit_height_std);
    fprintf(fileID, '  - Mean Summit Height (z_s):     %.4f\n', analysis_results.asperity.mean_summit_height);
end
fprintf(fileID, '==============================================================\n\n');
% --- Close the file if it was opened ---
if fileID ~= 1
    fclose(fileID);
    fprintf('Summary report saved to: %s\n', full_filepath);
end
end