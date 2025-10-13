# RSGA (Rough Surface Generator & Analyzer)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17341342.svg)](https://doi.org/10.5281/zenodo.17341342)


This is a modular and flexible MATLAB tool for generating and analyzing 2D randomly rough surfaces. It combines several methods into a single, configurable workflow.

## Features

- **Multiple Generation Methods**:
  - **`ACF`**: Generates a surface by filtering random noise. The spatial correlation is defined by an Autocovariance Function (ACF) model.
  - **`PSD_FFT`**: Generates a surface from a Power Spectral Density (PSD) definition using a fast, vectorized `ifft2` approach.
  - **`PSD_SUM`**: Generates a surface from a PSD using a slower, explicit summation method for direct comparison with older scripts.
- **Flexible Models**:
  - **ACF Models**: Includes Gaussian (`gaussian.m`) and Exponential (`exponential.m`) correlation models.
  - **PSD Models**: Includes a fractal model with a roll-off frequency (`simple_rolloff.m`) and the K-correlation model (`k_correlation.m`).
  - **Height Distributions**: Supports different initial noise profiles, such as standard Gaussian (`randn_noise.m`).
- **Comprehensive Analysis**:
  - **Statistical**: Calculates RMS Height (`Sq`), Skewness (`Ssk`), Kurtosis (`Sku`), and RMS Gradient (Sdq).
  - **Spectral**: Calculates spectral moments (`m0`, `m2`, `m4`) from both real-space gradients and by integrating the numerical PSD for self-consistency.
  - **Correlation**: Calculates the Autocovariance Function (ACF) and Correlation Lengths (`Clx`, `Cly`) for both isotropic and anisotropic surfaces.
  - **Hybrid parameter $\rho$** = (Sq * Sdq)/Cl [[Solhjoo and Vakis, Tribology International 115 (2017) 165-178]](https://pure.rug.nl/ws/portalfiles/portal/47462446/Chapter_5.pdf).
  - **Multi-Asperity**: Calculates Greenwood-Williamson style parameters derived from spectral moments, including summit density, mean summit radius, and more.
  - **PSD Comparison**: Includes two different radial averaging algorithms (`radialavg.m` and `RadiAve.m`) for comparison.
- **Save/Load Workflow**:
  - **Generate and Analyze**: Generate a new surface, analyze it, and save the complete session (config, data, results) to a `.mat` file in an `output` folder.
  - **Analyze Only**: Load a previously saved session file OR a simple `.mat` file containing only `x, y, f` data and run a full analysis.
- **Configurable**: All parameters and workflow options are controlled from a single, well-commented script, `main_runner.m`.

## File Structure

```
/SurfaceGeneratorTool/
|
|-- main_runner.m                 % The main script to configure and run.
|
|-- generate_surface.m            % Master function for surface generation.
|-- analyze_surface.m             % Master function for surface analysis.
|-- report_results.m              % Handles all plotting and command window output.
|
|-- /+psd_models/                 % Package for Power Spectral Density models.
|   |-- simple_rolloff.m
|   |-- k_correlation.m
|
|-- /+acf_models/                  % Package for Autocovariance Function models.
|   |-- gaussian.m
|   |-- exponential.m
|
|-- /+distributions/               % Package for initial noise distributions.
|   |-- randn_noise.m
|
|-- /+dist_pdf/                    % Package for analytical PDFs (for plotting).
|   |-- gaussian.m
|
|-- /utils/                        % Folder for helper functions.
|   |-- radialavg.m
|   |-- RadiAve.m
|   |-- solve_psd_constraint.m
|   |-- load_session.m
|   |-- save_session.m
|
|-- /output/                       % Default folder for saved .mat session files.
|
|-- README.md                      % This file.
```

## How to Use

All settings are controlled in the `main_runner.m` script.

### Workflow 1: Generate a New Surface

1.  Open `main_runner.m`.
2.  Set `config.workflow.mode = 'Generate and Analyze';`.
3.  Choose a generation method in `config.generation.method` (e.g., `'PSD_FFT'`, `'ACF'`).
4.  Adjust the grid settings and method-specific parameters as needed. For example, to generate a surface with an RMS height of 5.0 using the `simple_rolloff` PSD model:
    ```matlab
    config.generation.method = 'PSD_FFT';
    config.psd.model = @psd_models.simple_rolloff;
    config.psd.constraint_mode = 'm0';
    config.psd.target_rms_height = 5.0;
    config.psd.psd_slope = -2.5;
    config.psd.corr_lengths = 20.0;
    ```
5.  Set the desired analysis and plotting switches to `true` or `false`.
6.  Run the script.
7.  Plots will be displayed, a summary will be printed to the command window, and a `.mat` file with the complete session will be saved to the `output` folder.

### Workflow 2: Analyze an Existing Surface

This mode can load two types of files: a full session file created by this tool, or a simple `.mat` file from another source.

#### A) Analyzing a Full Session File

1.  Set `config.workflow.mode = 'Analyze Only';`.
2.  Set `config.workflow.loadFile` to the name of the file you want to load (e.g., `'output/GeneratedSurface_20251003_233000.mat'`).
3.  Set `config.workflow.overwrite_config_on_load = true;`. This will use the plotting and analysis switches from your current `main_runner.m` file, giving you control over the output.
4.  Run the script. The tool will skip generation, load the file, and run the requested analysis and reporting.

#### B) Analyzing a Simple `x, y, f` File

1.  Set `config.workflow.mode = 'Analyze Only';`.
2.  Set `config.workflow.loadFile` to the name of your simple `.mat` file (e.g., `'outSurf.mat'`).
3.  In the `config.workflow.loadSimpleFile.varNames` section, change the strings to match the variable names inside your `.mat` file. For `outSurf.mat`, the defaults are correct:
    ```matlab
    config.workflow.loadSimpleFile.varNames.x = 'x';
    config.workflow.loadSimpleFile.varNames.y = 'y';
    config.workflow.loadSimpleFile.varNames.f = 'f';
    ```
4.  Run the script. The tool will load the geometry, run a full analysis, and report the results. Since no target parameters exist, the report will not show comparisons.
