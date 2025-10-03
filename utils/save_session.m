function save_session(config, surface_data, analysis_results)
%SAVE_SESSION Saves the entire simulation session to a .mat file.
%   SAVE_SESSION(config, surface_data, analysis_results) bundles all
%   input and output data into a master struct and saves it.

disp('--- Saving Session ---');

% 1. Prepare a clean surface struct for saving, as requested.
save_surface_data = struct();
save_surface_data.x = surface_data.x_coords; % Already 1xN
save_surface_data.y = surface_data.y_coords; % Already 1xN
save_surface_data.f = surface_data.surface_h;  % Using 'f' for compatibility with user scripts
save_surface_data.generation_params = surface_data.generation_params;

% 2. Bundle all data into a single master struct.
session_data = struct();
session_data.config = config;
session_data.surface_data = save_surface_data;
session_data.analysis_results = analysis_results;

% 3. Check for/create the output directory.
output_folder = config.output.folder;
if ~exist(output_folder, 'dir')
    disp(['  - Output folder ''', output_folder, ''' not found. Creating it...']);
    mkdir(output_folder);
end

% 4. Generate a unique filename.
if isfield(config.output, 'filename_prefix') && ~isempty(config.output.filename_prefix)
    prefix = config.output.filename_prefix;
else
    prefix = 'SurfaceGenRun';
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = [prefix, '_', timestamp, '.mat'];
full_filepath = fullfile(output_folder, filename); % Create full path

% 5. Save the data.
try
    save(full_filepath, 'session_data');
    fprintf('Session successfully saved to: %s\n', full_filepath);
catch ME
    warning('Failed to save session data. Error: %s', ME.message);
end

end