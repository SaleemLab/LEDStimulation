function trials = crateTrialsStruct_flicker(nidq, csvFilePath)
%crateTrialsStruct_flicker Aligns trial data and calculates flicker phase vectors at 1 kHz.
%
%   This function intelligently aligns trials based on stimulus frequency, making it robust
%   to missing trials. It calculates the continuous phase vector in radians
%   for two color channels at a resampled rate of 1 kHz.
%
%   SYNTAX:
%   trials = crateTrialsStruct_flicker(nidq, csvFilePath)
%
%   INPUTS:
%   nidq         - A struct from extractNidqData, containing 'stimON', 'stimCHANGE', etc.
%   csvFilePath  - Path to the CSV file with trial parameters. May have a header.
%
%   OUTPUT:
%   trials       - A struct array where each element represents an aligned trial,
%                  containing parameters, measured data, and phase vectors.

%% 1. Load Data Sources
% Use imecTime if available, otherwise default to NIDQ time
if isfield(nidq, 'imecTime') && ~isempty(nidq.imecTime)
    timeVector = nidq.imecTime;
    timeBaseSource = 'imecTime';
else
    timeVector = nidq.time;
    timeBaseSource = 'NIDQ time';
end

% Load trial parameters from CSV, automatically handling a potential header
if ~isfile(csvFilePath)
    error('CSV file not found at: %s', csvFilePath);
end
opts = detectImportOptions(csvFilePath);
% Standardize variable names to use 'green' and 'uv'
opts.VariableNames = {'stimType', 'duration_ms', 'frequency_param', 'phase_green', 'phase_uv', ...
    'contrast_green', 'contrast_uv', 'bonsaiTime', 'arduinoTime'};
trialParams = readtable(csvFilePath, opts);
num_csv_trials = height(trialParams);
fprintf('Loaded %d trials from parameter file.\n', num_csv_trials);

%% 2. Pre-process NIDQ Data to Find All Potential Trials
% Find trial boundaries from the 'stimON' signal
on_indices  = find(diff(nidq.stimON) == 1) + 1;
off_indices = find(diff(nidq.stimON) == -1) + 1;

% Clean up incomplete trials at the edges
if ~isempty(off_indices) && ~isempty(on_indices)
    if off_indices(1) < on_indices(1), off_indices(1) = []; end
    if ~isempty(on_indices) && on_indices(end) > off_indices(end), on_indices(end) = []; end
end

num_nidq_trials = length(on_indices);
fprintf('Found %d potential trials in NIDQ data.\n', num_nidq_trials);

% Create a temporary structure with measured data for each NIDQ trial
nidq_trials_info(num_nidq_trials,1) = struct('startTime',[],'endTime',[],'measured_frequency',[]);
for i = 1:num_nidq_trials
    nidq_trials_info(i).startTime = timeVector(on_indices(i));
    nidq_trials_info(i).endTime = timeVector(off_indices(i));
    
    % Extract stimCHANGE snippet and calculate measured frequency
    time_mask = (timeVector >= nidq_trials_info(i).startTime & timeVector < nidq_trials_info(i).endTime);
    stimCHANGE_snippet = nidq.stimCHANGE(time_mask);
    toggle_indices = find(diff(stimCHANGE_snippet) ~= 0) + 1;
    
    if length(toggle_indices) < 2
        nidq_trials_info(i).measured_frequency = NaN;
    else
        toggle_times = timeVector(on_indices(i) + toggle_indices - 1);
        nidq_trials_info(i).measured_frequency = 1 / mean(diff(toggle_times));
    end
end

%% 3. Align NIDQ Trials with CSV Parameters using Frequency
matched_trials = struct(); % Initialize an empty struct for matched trials
nidq_idx = 1;
csv_idx = 1;
freq_tolerance = 1.0; % Allow 1 Hz difference for a match
trial_count = 0;

while nidq_idx <= num_nidq_trials && csv_idx <= num_csv_trials
    
    % --- Core Alignment Logic ---
    freq_nidq = nidq_trials_info(nidq_idx).measured_frequency;
    freq_csv = trialParams.frequency_param(csv_idx);
    is_match_flag = ~isnan(freq_nidq) && abs(freq_nidq - freq_csv) < freq_tolerance;

    if is_match_flag
        % SUCCESSFUL MATCH: Frequencies align. Append to our list.
        trial_count = trial_count + 1;
        matched_trials = append_matched_trial(matched_trials, trial_count, nidq_trials_info(nidq_idx), trialParams(csv_idx, :), timeBaseSource);
        nidq_idx = nidq_idx + 1;
        csv_idx = csv_idx + 1;
    else
        % MISMATCH: Look ahead to determine which trial is extra.
        lookahead_nidq_matches = false;
        if nidq_idx + 1 <= num_nidq_trials
            next_freq_nidq = nidq_trials_info(nidq_idx + 1).measured_frequency;
            lookahead_nidq_matches = ~isnan(next_freq_nidq) && abs(next_freq_nidq - freq_csv) < freq_tolerance;
        end

        if lookahead_nidq_matches
            % The current NIDQ trial is extra/unmatched. Skip it.
            nidq_idx = nidq_idx + 1;
        else
            % Assume the CSV trial is the extra one. Skip it.
            csv_idx = csv_idx + 1;
        end
    end
end
trials = matched_trials; % Assign the final struct
fprintf('\n--- Alignment Complete: %d trials successfully matched. ---\n', length(trials));
end

%% --- Helper function for matched trials ---
function s = append_matched_trial(s, idx, nidq_info, csv_row, timeBaseSource)
    s(idx).timeBase = timeBaseSource;
    s(idx).startTime = nidq_info.startTime;
    s(idx).endTime = nidq_info.endTime;
    s(idx).duration_measured_s = s(idx).endTime - s(idx).startTime;
    s(idx).duration_param_s = csv_row.duration_ms / 1000;
    s(idx).frequency_measured = nidq_info.measured_frequency;
    s(idx).frequency_param = csv_row.frequency_param;
    s(idx).stimType = csv_row.stimType{1};
    s(idx).contrast_green = csv_row.contrast_green;
    s(idx).contrast_uv = csv_row.contrast_uv;
    s(idx).phase_green_param_deg = csv_row.phase_green;
    s(idx).phase_uv_param_deg = csv_row.phase_uv;
    
    % --- Phase Vector Calculation at 1 kHz ---
    if isnan(s(idx).frequency_measured)
        s(idx).time_1khz = [];
        s(idx).phase_green_rad = [];
        s(idx).phase_uv_rad = [];
    else
        resampleFs = 1000;
        dt = 1/resampleFs;
        s(idx).time_1khz = (s(idx).startTime : dt : s(idx).endTime)';
        time_elapsed = s(idx).time_1khz - s(idx).startTime;
        unwrapped_phase_rad = 2 * pi * s(idx).frequency_measured * time_elapsed;
        phase_green_offset_rad = deg2rad(s(idx).phase_green_param_deg);
        phase_uv_offset_rad = deg2rad(s(idx).phase_uv_param_deg);
        s(idx).phase_green_rad = mod(unwrapped_phase_rad + phase_green_offset_rad, 2*pi);
        s(idx).phase_uv_rad    = mod(unwrapped_phase_rad + phase_uv_offset_rad, 2*pi);
    end

    % Validation
    if abs(s(idx).duration_measured_s - s(idx).duration_param_s) > 0.05
        warning('Trial %d (Matched): Duration mismatch. Measured: %.3fs, Param: %.3fs', idx, s(idx).duration_measured_s, s(idx).duration_param_s);
    end
end