function surface_data = generate_surface(config)
%GENERATE_SURFACE Generates a 2D rough surface based on the specified method.
%   surface_data = GENERATE_SURFACE(config) takes a configuration struct
%   and returns a struct containing the generated surface and its grid.
%
%   The function acts as a dispatcher, calling the appropriate sub-method
%   based on config.generation.method ('PSD' or 'ACF').
%
%   Output surface_data struct contains:
%       .x_coords       (1xNx vector of x-coordinates)
%       .y_coords       (1xNy vector of y-coordinates)
%       .X_grid         (NyxNx matrix of x-coordinates)
%       .Y_grid         (NyxNx matrix of y-coordinates)
%       .surface_h      (NyxNx matrix of surface heights)

disp('--- Starting Surface Generation ---');
validate_config(config);
num_points = config.grid.num_points;
point_spacing = config.grid.point_spacing;
Nx = num_points(1); Ny = num_points(2);
% dx = point_spacing(1); dy = point_spacing(2);

% --- Initialize the output struct ONCE at the beginning ---
lateral_length = (num_points - 1) .* point_spacing;
x_coords = linspace(-lateral_length(1)/2, lateral_length(1)/2, Nx);
y_coords = linspace(-lateral_length(2)/2, lateral_length(2)/2, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

surface_data = struct('x_coords', x_coords, 'y_coords', y_coords, ...
    'X_grid', X_grid, 'Y_grid', Y_grid);

switch config.generation.method
    case 'ACF'
        % This case is correct and remains unchanged
        disp("Using 'ACF' method...");
        if config.grid.is_anisotropic && isscalar(config.acf.corr_lengths)
            warning('Configuration specifies an ANISOTROPIC grid, but a single value was given for corr_lengths. Proceeding with ISOTROPIC correlation.');
        end
        raw_noise = config.acf.initial_noise_dist(num_points);
        acf_filter = config.acf.model(X_grid, Y_grid, config.acf.corr_lengths);
        correlated_surface = real(ifft2(fft2(raw_noise) .* fft2(acf_filter)));
        correlated_surface = correlated_surface - mean(correlated_surface(:));
        current_rms = std(correlated_surface(:));
        scaling_factor = config.acf.target_rms_height / current_rms;
        surface_h = correlated_surface * scaling_factor;
        disp(['  - Surface generated with RMS height: ', num2str(std(surface_h(:))), ' units.']);

    case 'PSD_FFT'
        disp("Using 'PSD' method (Fast Fourier Transform)...");

        Lx = lateral_length(1); Ly = lateral_length(2);
        qx_vec = 2*pi * (-Nx/2 : Nx/2-1) / Lx;
        qy_vec = 2*pi * (-Ny/2 : Ny/2-1) / Ly;
        [Qx_grid, Qy_grid] = meshgrid(qx_vec, qy_vec);

        % --- Handle different model types ---
        model_name = func2str(config.psd.model);
        if contains(model_name, 'k_correlation')
            % This model calculates its own amplitude.
            disp(['  - Generating 2D PSD with ''', model_name, ''' model.']);
            psd_2d = config.psd.model(Qx_grid, Qy_grid, [], config.psd);
            
            % For reporting, create dummy generation_params
            surface_data.generation_params.q_lower = min(abs(qx_vec(qx_vec~=0)));
            surface_data.generation_params.q_upper = max(abs(qx_vec));
            surface_data.generation_params.q_rolloff = NaN;
            surface_data.generation_params.C0 = NaN;
        elseif contains(model_name, 'simple_rolloff')
            % For other models, solve for the C0 amplitude.
            dqx = 2*pi / Lx;
            dqy = 2*pi / Ly;
            L_eff = max(lateral_length);
            delta_eff = min(point_spacing);
            q_lower = 2*pi / L_eff;
            q_upper = pi / delta_eff;
            q_grids = struct('Qx_grid', Qx_grid, 'Qy_grid', Qy_grid, 'dqx', dqx, 'dqy', dqy, 'q_lower', q_lower, 'q_upper', q_upper);
            
            disp(['  - Solving for PSD parameters with constraint mode: ''', config.psd.constraint_mode, '''...']);
            [C0, solved_params] = solve_psd_constraint(config, q_grids);
            
            disp(['  - Generating 2D PSD with ''', func2str(config.psd.model), ''' model.']);
            psd_2d = config.psd.model(Qx_grid, Qy_grid, C0, solved_params);
            
            surface_data.generation_params.q_lower = solved_params.q_lower;
            surface_data.generation_params.q_upper = solved_params.q_upper;
            surface_data.generation_params.q_rolloff = solved_params.q_rolloff;
            surface_data.generation_params.C0 = C0;
        else
            error('Unsupported PSD model: "%s". The generation logic in generate_surface.m needs to be updated to handle this model.', model_name);
        end
        
        % Set amplitude proportional to sqrt(PSD) to get the correct relative power (shape)
        amplitude_2d = sqrt(psd_2d);

        if strcmp(config.psd.amplitude_method, 'fft_filter')
            disp('  - Applying stochastic amplitudes (FFT Filtering Method).');
            amplitude_2d = amplitude_2d .* randn(Ny, Nx);
        end

        random_phases = 2 * pi * rand(Ny, Nx);
        conj_symm_phases = random_phases - rot90(random_phases, 2);

        fourier_space_surface = amplitude_2d .* exp(1i * conj_symm_phases);

        % disp('  - Performing Inverse Fast Fourier Transform...');
        surface_h_unscaled = real(ifft2(ifftshift(fourier_space_surface)));

        % Normalize the generated surface and scale it to the exact target RMS height
        current_rms = std(surface_h_unscaled(:));
        target_rms = config.psd.target_rms_height;

        % Avoid division by zero if the surface is flat
        if current_rms > 0
            scaling_factor = target_rms / current_rms;
            surface_h = surface_h_unscaled * scaling_factor;
        else
            surface_h = surface_h_unscaled; % Already zero
        end

        % surface_data.generation_params.q_lower = q_lower;
        % surface_data.generation_params.q_upper = q_upper;
        % surface_data.generation_params.q_rolloff = solved_params.q_rolloff;
        % surface_data.generation_params.C0 = C0;


    case 'PSD_SUM'
        disp("Using 'PSD' method (Explicit Sums)...");

        L_eff = max(lateral_length);
        delta_eff = min(point_spacing);
        q_lower = 2*pi / L_eff;
        q_upper = pi / delta_eff;
        grid_params = struct('q_lower', q_lower, 'q_upper', q_upper);

        disp(['  - Solving for PSD amplitude with constraint mode: ''', config.psd.constraint_mode, '''...']);
        [C0, solved_params] = solve_psd_constraint(config, grid_params);
        q_rolloff = solved_params.q_rolloff;
        psd_slope = solved_params.psd_slope;

        q_step = q_lower;
        qx_vec = (ceil(-q_upper/q_step):1:floor(q_upper/q_step)) * q_step;
        qy_vec = qx_vec;
        Nqx = length(qx_vec);
        Nqy = length(qy_vec);

        phi_u = rand(floor(Nqy/2), Nqx);
        phi_d = -rot90(phi_u, 2);
        if mod(Nqy, 2) == 0
            random_phases = [phi_u; phi_d];
        else
            phi_c_rand = rand(1, floor(Nqx/2));
            phi_c_symm = -rot90(phi_c_rand, 2);
            phi_c = [phi_c_symm, 0, phi_c_rand];
            random_phases = [phi_u; phi_c; phi_d];
        end
        random_phases = 2 * pi * random_phases;

        disp('  - Generating surface via explicit summation loop (this may be slow)...');
        surface_h = zeros(Ny, Nx, 'like', 1i);
        % amplitude_sqrt_C0 = sqrt(C0);

        % for i = 1:Nqx
        %     for j = 1:Nqy
        %         q_mag = sqrt(qx_vec(i)^2 + qy_vec(j)^2);
        %         amplitude = 0;
        % 
        %         if (q_mag >= q_lower) && (q_mag < q_rolloff)
        %             amplitude = amplitude_sqrt_C0;
        %         elseif (q_mag >= q_rolloff) && (q_mag <= q_upper)
        %             amplitude = amplitude_sqrt_C0 * (q_mag / q_rolloff)^(psd_slope/2);
        %         end
        % 
        %         if amplitude > 0
        %             phase = random_phases(j,i);
        %             surface_h = surface_h + amplitude * exp(1i * (qx_vec(i) * X_grid + qy_vec(j) * Y_grid + phase));
        %         end
        %     end
        % end
        model_name = func2str(config.psd.model);

        for i = 1:Nqx
            for j = 1:Nqy
                q_mag = sqrt(qx_vec(i)^2 + qy_vec(j)^2);
                amplitude = 0;

                if contains(model_name, 'simple_rolloff')
                    if (q_mag >= q_lower) && (q_mag < q_rolloff)
                        amplitude = sqrt(C0);
                    elseif (q_mag >= q_rolloff) && (q_mag <= q_upper)
                        amplitude = sqrt(C0 * (q_mag / q_rolloff)^psd_slope);
                    end

                elseif contains(model_name, 'k_correlation')
                    h0 = config.psd.target_rms_height;
                    cl = config.psd.corr_lengths;
                    m = config.psd.psd_slope;

                    term = 1 - (cl^2 / (m + 2)) * q_mag^2;
                    if term > 0
                        psd_val = ((cl*h0)^2 / (2*pi)) * term^(m/2);
                        amplitude = sqrt(psd_val);
                    end
                end

                if amplitude > 0
                    phase = random_phases(j,i);
                    surface_h = surface_h + amplitude * exp(1i * (qx_vec(i) * X_grid + qy_vec(j) * Y_grid + phase));
                end
            end
        end
        surface_h = real(surface_h * q_lower);

        disp(['  - Surface generated with RMS height: ', num2str(std(surface_h(:))), ' units.']);

        surface_data.generation_params.q_lower = q_lower;
        surface_data.generation_params.q_upper = q_upper;
        surface_data.generation_params.q_rolloff = solved_params.q_rolloff;
        surface_data.generation_params.C0 = C0;

    % case 'Iterative'
    %     error("Generation method 'Iterative' is not yet implemented.");
end

% surface_data.x_coords = x_coords;
% surface_data.y_coords = y_coords;
% surface_data.X_grid = X_grid;
% surface_data.Y_grid = Y_grid;
surface_data.surface_h = surface_h;

disp('--- Surface Generation Complete ---');
end

function validate_config(config)
% Helper function for input validation
assert(isfield(config, 'grid'), 'Config must contain a "grid" field.');
assert(isfield(config.grid, 'num_points') && isvector(config.grid.num_points) && numel(config.grid.num_points) == 2, ...
    'config.grid.num_points must be a 1x2 vector.');
assert(isfield(config.grid, 'point_spacing') && isvector(config.grid.point_spacing) && numel(config.grid.point_spacing) == 2, ...
    'config.grid.point_spacing must be a 1x2 vector.');
assert(isfield(config.grid, 'is_anisotropic') && islogical(config.grid.is_anisotropic), ...
    'config.grid.is_anisotropic must be true or false.');

assert(isfield(config, 'generation') && isfield(config.generation, 'method'), ...
    'Config must specify a generation method in config.generation.method.');

method = config.generation.method;
valid_methods = {'ACF', 'PSD_FFT', 'PSD_SUM'};
assert(ismember(method, valid_methods), 'config.generation.method must be one of: ''ACF'', ''PSD_FFT'', ''PSD_SUM''.');

% --- Method-specific checks ---
if strcmp(method, 'ACF')
    assert(isfield(config, 'acf'), 'For ACF method, config must contain an "acf" field.');
    assert(isfield(config.acf, 'corr_lengths'), 'config.acf must contain "corr_lengths".');

    if config.grid.is_anisotropic
        % If anisotropic is intended, check if input is a vector.
        % If it's a scalar, we will allow it but warn the user later.
        % If it's neither a scalar nor a 1x2 vector, it's an error.
        if ~isscalar(config.acf.corr_lengths)
            assert(isvector(config.acf.corr_lengths) && numel(config.acf.corr_lengths) == 2, ...
                'For ANISOTROPIC generation, config.acf.corr_lengths must be a 1x2 vector [cl_x, cl_y].');
        end
    else
        % If isotropic is intended, it MUST be a scalar.
        assert(isscalar(config.acf.corr_lengths), ...
            'For ISOTROPIC generation, config.acf.corr_lengths must be a single scalar value.');
    end
end

% (Add checks for 'PSD' and 'Iterative' methods here as they are implemented)
end