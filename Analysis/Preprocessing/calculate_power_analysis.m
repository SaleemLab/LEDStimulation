function [units_power_analysis, f] = calculate_power_analysis(units_struct, allIntervalData, uniqueFreqs, chronux_params, target_duration_ms, duration_tolerance_ms)
% CALCULATE_POWER_ANALYSIS Performs batch power spectrum analysis for all units.
%
%   ... (function description) ...
%
%   INPUTS:
%   - units_struct:     ...
%   - allIntervalData:  ...
%   - uniqueFreqs:      ...
%   - chronux_params:   ...
%   - target_duration_ms: The target trial duration in ms (e.g., 2000).
%   - duration_tolerance_ms: Tolerance for filtering trials (e.g., 50).
%
%   OUTPUTS:
%   - units_power_analysis: ...
%   - f:                    ...
%

% --- 1. Initialization ---
nUnits = numel(units_struct);
nLums = size(allIntervalData, 2);
nContrasts = size(allIntervalData, 1);
nFreqs = numel(uniqueFreqs) - 1; % -1 because baseline (200Hz) is separate

f = []; % Frequency vector, to be populated once

% --- VALIDATE TAPER PARAMS ---
if ~isfield(chronux_params, 'tapers') || ...
   ~isnumeric(chronux_params.tapers) || ...
   numel(chronux_params.tapers) ~= 2
    error('calculate_power_analysis:InvalidTapers', ...
          'Input chronux_params.tapers must be in the [TW, K] format (e.g., [2 4]).');
end
% --- END VALIDATION ---

% --- VALIDATE DURATION PARAMS ---
if nargin < 5
    error('calculate_power_analysis:MissingInput', 'target_duration_ms (in ms) is a required input.');
end
if nargin < 6
    duration_tolerance_ms = 50; % Default tolerance of 50ms
    fprintf('Using default duration tolerance of 50ms.\n');
end
% --- END VALIDATION ---


% Store the original tapers, which we will use to calculate W and p
default_tapers = chronux_params.tapers; % e.g., [2 4]
TW = default_tapers(1);
K = default_tapers(2);
p = 2*TW - K;

% --- Create a FIXED taper struct for all calls ---
% We will pass the [TW, K] tapers directly, and also pass a 't'
% field, which explicitly defines the time grid. This is the
% most robust way to ensure all calls to mtspectrumpt use the
% same FFT length, even for empty trials.
T_sec = target_duration_ms / 1000.0;
fixed_taper_params = chronux_params;
% Ensure tapers are in [TW, K] format (which they are, via validation)
fixed_taper_params.tapers = default_tapers; 
% Add the explicit time grid
fixed_taper_params.t = 0:(1/fixed_taper_params.Fs):T_sec;
% --- End fixed taper creation ---


% Pre-allocate the output struct array
units_power_analysis(nUnits) = struct();

for iunit = 1:nUnits
    units_power_analysis(iunit).unit_id = units_struct(iunit).cluster_id;
    % Pre-allocate the .powerAnalysis struct and its fields
    units_power_analysis(iunit).powerAnalysis = struct();
    units_power_analysis(iunit).powerAnalysis.baselinePower_stat = cell(1, nLums);
    units_power_analysis(iunit).powerAnalysis.baselinePower_run  = cell(1, nLums);
    units_power_analysis(iunit).powerAnalysis.freq_stat = [];
    units_power_analysis(iunit).powerAnalysis.freq_run  = [];
    
    % Pre-allocate all condition cells
    units_power_analysis(iunit).powerAnalysis.power_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.power_baselineNorm_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1norm_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1_double_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1norm_double_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1_half_stat = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1norm_half_stat = cell(nContrasts, nLums, nFreqs);

    units_power_analysis(iunit).powerAnalysis.power_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.power_baselineNorm_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1norm_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1_double_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1norm_double_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1_half_run = cell(nContrasts, nLums, nFreqs);
    units_power_analysis(iunit).powerAnalysis.powerF1norm_half_run = cell(nContrasts, nLums, nFreqs);
end

fprintf('Starting power analysis for %d units...\n', nUnits);

% --- 2. Get baseline data (pre-loop) ---
stat_baseline_alllumi = cell(1, nLums);
run_baseline_alllumi = cell(1, nLums);
for iLum = 1:nLums
    stat_baseline_alllumi{iLum} = cat(1, allIntervalData{1, iLum}.stat_baseline, allIntervalData{2, iLum}.stat_baseline);
    run_baseline_alllumi{iLum} = cat(1, allIntervalData{1, iLum}.run_baseline, allIntervalData{2, iLum}.run_baseline);
end

% --- 3. Main Unit Loop ---
for iunit = 1:nUnits
    if mod(iunit, 10) == 0
        fprintf('  Processing unit %d of %d\n', iunit, nUnits);
    end
    
    unit_spike_times = units_struct(iunit).spike_times .* 1000; 

    for iLum = 1:nLums
        % --- Stat Baseline ---
        intervals_stat = stat_baseline_alllumi{iLum};
        
        % --- FILTER INTERVALS BY DURATION ---
        if ~isempty(intervals_stat)
            durations_ms = intervals_stat(:, 2) - intervals_stat(:, 1);
            good_idx = abs(durations_ms - target_duration_ms) <= duration_tolerance_ms;
            intervals_stat = intervals_stat(good_idx, :); % Keep only good intervals
        end
        % --- END FILTER ---
        
        [S_stat, f_stat, ~, ~] = calculate_spectrum_from_intervals(unit_spike_times, intervals_stat, fixed_taper_params);
        units_power_analysis(iunit).powerAnalysis.baselinePower_stat{iLum} = S_stat;
        
        % --- Run Baseline ---
        intervals_run = run_baseline_alllumi{iLum};
        
        % --- FILTER INTERVALS BY DURATION ---
        if ~isempty(intervals_run)
            durations_ms = intervals_run(:, 2) - intervals_run(:, 1);
            good_idx = abs(durations_ms - target_duration_ms) <= duration_tolerance_ms;
            intervals_run = intervals_run(good_idx, :);
        end
        % --- END FILTER ---
        
        [S_run, f_run, ~, ~] = calculate_spectrum_from_intervals(unit_spike_times, intervals_run, fixed_taper_params);
        units_power_analysis(iunit).powerAnalysis.baselinePower_run{iLum} = S_run;

        % Get the mean baseline power for normalization
        mean_baseline_stat = mean(units_power_analysis(iunit).powerAnalysis.baselinePower_stat{iLum}, 2);
        mean_baseline_run = mean(units_power_analysis(iunit).powerAnalysis.baselinePower_run{iLum}, 2);

        mean_baseline_stat(mean_baseline_stat == 0) = NaN; % Prevents Inf
        mean_baseline_run(mean_baseline_run == 0) = NaN; % Prevents Inf

        % Store frequency vector (only need to do this once)
        if isempty(f)
            f = f_stat; 
            units_power_analysis(iunit).powerAnalysis.freq_stat = f_stat; % Store once
            units_power_analysis(iunit).powerAnalysis.freq_run  = f_run;  % Store once
        end
        
        for iContrast = 1:nContrasts
            for ifreq = 1:nFreqs
               
                
                % --- Stat Power Spectrum --- %
                intervals_stat = allIntervalData{iContrast, iLum}.stat_intervals{ifreq};
                
                % --- FILTER INTERVALS BY DURATION ---
                if ~isempty(intervals_stat)
                    durations_ms = intervals_stat(:, 2) - intervals_stat(:, 1);
                    good_idx = abs(durations_ms - target_duration_ms) <= duration_tolerance_ms;
                    intervals_stat = intervals_stat(good_idx, :);
                end
                % --- END FILTER ---
                
                [S_stat, ~, ~, ~] = calculate_spectrum_from_intervals(unit_spike_times, intervals_stat, fixed_taper_params);
                units_power_analysis(iunit).powerAnalysis.power_stat{iContrast, iLum, ifreq} = S_stat;
                
                % Normalize
                if ~isempty(S_stat) && size(S_stat, 2) > 0 && ~isempty(mean_baseline_stat)
                    S_norm_stat = S_stat ./ mean_baseline_stat;
                    units_power_analysis(iunit).powerAnalysis.power_baselineNorm_stat{iContrast, iLum, ifreq} = S_norm_stat;
                else
                    S_norm_stat = []; 
                    units_power_analysis(iunit).powerAnalysis.power_baselineNorm_stat{iContrast, iLum, ifreq} = [];
                end

                % Power at F1, F1/2 and 2*F1
                if ~isempty(S_stat) && size(S_stat, 2) > 0 && ~isempty(f)
                    target = uniqueFreqs(ifreq); %F1
                    [~, idx] = min(abs(f - target));
                    units_power_analysis(iunit).powerAnalysis.powerF1_stat{iContrast, iLum, ifreq} = S_stat(idx, :);
                    units_power_analysis(iunit).powerAnalysis.powerF1norm_stat{iContrast, iLum, ifreq} = S_norm_stat(idx, :);
                    
                    target  = uniqueFreqs(ifreq) * 2; %F1*2
                    [~, idx] = min(abs(f - target));
                    units_power_analysis(iunit).powerAnalysis.powerF1_double_stat{iContrast, iLum, ifreq} = S_stat(idx, :);
                    units_power_analysis(iunit).powerAnalysis.powerF1norm_double_stat{iContrast, iLum, ifreq} = S_norm_stat(idx, :);
                    
                    target  = uniqueFreqs(ifreq) / 2; %F1/2
                    [~, idx] = min(abs(f - target));
                    units_power_analysis(iunit).powerAnalysis.powerF1_half_stat{iContrast, iLum, ifreq} = S_stat(idx, :);
                    units_power_analysis(iunit).powerAnalysis.powerF1norm_half_stat{iContrast, iLum, ifreq} = S_norm_stat(idx, :);
                end
                
                % --- Run Power Spectrum --- %
                intervals_run = allIntervalData{iContrast, iLum}.run_intervals{ifreq};
                
                % --- FILTER INTERVALS BY DURATION ---
                if ~isempty(intervals_run)
                    durations_ms = intervals_run(:, 2) - intervals_run(:, 1);
                    good_idx = abs(durations_ms - target_duration_ms) <= duration_tolerance_ms;
                    intervals_run = intervals_run(good_idx, :);
                end
                % --- END FILTER ---
                
                [S_run, ~, ~, ~] = calculate_spectrum_from_intervals(unit_spike_times, intervals_run, fixed_taper_params);
                units_power_analysis(iunit).powerAnalysis.power_run{iContrast, iLum, ifreq} = S_run;

                % Normalize
                if ~isempty(S_run) && size(S_run, 2) > 0 && ~isempty(mean_baseline_run)
                    S_norm_run = S_run ./ mean_baseline_run;
                    units_power_analysis(iunit).powerAnalysis.power_baselineNorm_run{iContrast, iLum, ifreq} = S_norm_run;
                else
                    S_norm_run = [];
                    units_power_analysis(iunit).powerAnalysis.power_baselineNorm_run{iContrast, iLum, ifreq} = [];
                end
                
                % Power at F1, F1/2 and 2*F1
                if ~isempty(S_run) && size(S_run, 2) > 0 && ~isempty(f)
                    target = uniqueFreqs(ifreq); %F1
                    [~, idx] = min(abs(f - target));
                    units_power_analysis(iunit).powerAnalysis.powerF1_run{iContrast, iLum, ifreq} = S_run(idx, :);
                    units_power_analysis(iunit).powerAnalysis.powerF1norm_run{iContrast, iLum, ifreq} = S_norm_run(idx, :);
                    
                    target  = uniqueFreqs(ifreq) * 2; %F1*2
                    [~, idx] = min(abs(f - target));
                    units_power_analysis(iunit).powerAnalysis.powerF1_double_run{iContrast, iLum, ifreq} = S_run(idx, :);
                    units_power_analysis(iunit).powerAnalysis.powerF1norm_double_run{iContrast, iLum, ifreq} = S_norm_run(idx, :);
                    
                    target  = uniqueFreqs(ifreq) / 2; %F1/2
                    [~, idx] = min(abs(f - target));
                    units_power_analysis(iunit).powerAnalysis.powerF1_half_run{iContrast, iLum, ifreq} = S_run(idx, :);
                    units_power_analysis(iunit).powerAnalysis.powerF1norm_half_run{iContrast, iLum, ifreq} = S_norm_run(idx, :);
                end
            end
        end
    end
end

fprintf('Power analysis complete.\n');
end

