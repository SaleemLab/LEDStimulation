function trials = crateTrialsStruct_flicker(nidq, csvFilePath)
%crateTrialsStruct_flicker Aligns and merges trial data from NIDQ signals and a parameter file.
%
%   This function intelligently aligns trials based on stimulus frequency,
%   making it robust to missing trials in either the NIDQ data or the CSV
%   log. It automatically detects and skips a header row in the CSV file if
%   one is present.
%
%   SYNTAX:
%   trials = crateTrialsStruct_flicker(nidq, csvFilePath)
%
%   INPUTS:
%   nidq         - A struct from extractNidqData, containing 'stimON', 'stimCHANGE', etc.
%   csvFilePath  - Path to the CSV file with trial parameters. May have a header.
%
%   OUTPUT:
%   trials       - A struct array where each element is a trial. Includes a
%                  'match_status' field ('Success', 'Unmatched NIDQ', 'Unmatched CSV').

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
opts.VariableNames = {'stimType', 'duration_ms', 'frequency_param', 'phase1', 'phase2', ...
    'contrast1', 'contrast2', 'bonsaiTime', 'arduinoTime'};
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
trials = struct(); % Initialize the final output struct
trial_count = 0;
nidq_idx = 1;
csv_idx = 1;
freq_tolerance = 1.0; % Allow 1 Hz difference for a match

while nidq_idx <= num_nidq_trials || csv_idx <= num_csv_trials
    trial_count = trial_count + 1;
    
    % Check if we have run out of trials in one of the sources
    if nidq_idx > num_nidq_trials
        % Ran out of NIDQ trials; remaining CSV trials are unmatched
        trials = append_unmatched_csv_trial(trials, trial_count, trialParams(csv_idx, :), timeBaseSource);
        csv_idx = csv_idx + 1;
        continue;
    end
    if csv_idx > num_csv_trials
        % Ran out of CSV trials; remaining NIDQ trials are unmatched
        trials = append_unmatched_nidq_trial(trials, trial_count, nidq_trials_info(nidq_idx), nidq, timeVector, timeBaseSource);
        nidq_idx = nidq_idx + 1;
        continue;
    end

    % --- Core Alignment Logic ---
    freq_nidq = nidq_trials_info(nidq_idx).measured_frequency;
    freq_csv = trialParams.frequency_param(csv_idx);
    
    is_match = ~isnan(freq_nidq) && abs(freq_nidq - freq_csv) < freq_tolerance;

    if is_match
        % SUCCESSFUL MATCH: Frequencies align.
        trials = append_matched_trial(trials, trial_count, nidq_trials_info(nidq_idx), trialParams(csv_idx, :), nidq, timeVector, timeBaseSource);
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
            % The current NIDQ trial is extra/unmatched.
            trials = append_unmatched_nidq_trial(trials, trial_count, nidq_trials_info(nidq_idx), nidq, timeVector, timeBaseSource);
            nidq_idx = nidq_idx + 1; % Skip this NIDQ trial and try again
        else
            % Assume the CSV trial is the extra one.
            trials = append_unmatched_csv_trial(trials, trial_count, trialParams(csv_idx, :), timeBaseSource);
            csv_idx = csv_idx + 1; % Skip this CSV trial and try again
        end
    end
end

% Final report
matched_count = sum(strcmp({trials.match_status}, 'Success'));
unmatched_nidq_count = sum(strcmp({trials.match_status}, 'Unmatched NIDQ'));
unmatched_csv_count = sum(strcmp({trials.match_status}, 'Unmatched CSV'));

fprintf(['\n--- Alignment Complete ---\n' ...
    'Successfully matched: %d trials\n' ...
    'Unmatched NIDQ trials: %d\n' ...
    'Unmatched CSV trials: %d\n' ...
    'Total entries in struct: %d\n'], ...
    matched_count, unmatched_nidq_count, unmatched_csv_count, length(trials));
end

% --- Helper functions to keep the main loop clean ---
function s = append_matched_trial(s, idx, nidq_info, csv_row, nidq, timeVector, timeBaseSource)
    s(idx).trial_idx = idx;
    s(idx).match_status = 'Success';
    s(idx).timeBase = timeBaseSource;
    s(idx).startTime = nidq_info.startTime;
    s(idx).endTime = nidq_info.endTime;
    s(idx).duration_measured_s = s(idx).endTime - s(idx).startTime;
    s(idx).duration_param_s = csv_row.duration_ms / 1000;
    s(idx).frequency_measured = nidq_info.measured_frequency;
    s(idx).frequency_param = csv_row.frequency_param;
    s(idx).stimType = csv_row.stimType{1};
    s(idx).contrast1 = csv_row.contrast1;
    s(idx).contrast2 = csv_row.contrast2;
    
    % Add data snippets
    time_mask = (timeVector >= s(idx).startTime & timeVector < s(idx).endTime);
    s(idx).time = timeVector(time_mask);
    s(idx).stimCHANGE = nidq.stimCHANGE(time_mask);
    
    % Validation
    if abs(s(idx).duration_measured_s - s(idx).duration_param_s) > 0.05
        warning('Trial %d (Matched): Duration mismatch. Measured: %.3fs, Param: %.3fs', idx, s(idx).duration_measured_s, s(idx).duration_param_s);
    end
end

function s = append_unmatched_nidq_trial(s, idx, nidq_info, nidq, timeVector, timeBaseSource)
    s(idx).trial_idx = idx;
    s(idx).match_status = 'Unmatched NIDQ';
    s(idx).timeBase = timeBaseSource;
    s(idx).startTime = nidq_info.startTime;
    s(idx).endTime = nidq_info.endTime;
    s(idx).duration_measured_s = s(idx).endTime - s(idx).startTime;
    s(idx).duration_param_s = NaN;
    s(idx).frequency_measured = nidq_info.measured_frequency;
    s(idx).frequency_param = NaN;
    s(idx).stimType = 'N/A';
    s(idx).contrast1 = NaN;
    s(idx).contrast2 = NaN;
    time_mask = (timeVector >= s(idx).startTime & timeVector < s(idx).endTime);
    s(idx).time = timeVector(time_mask);
    s(idx).stimCHANGE = nidq.stimCHANGE(time_mask);
end

function s = append_unmatched_csv_trial(s, idx, csv_row, timeBaseSource)
    s(idx).trial_idx = idx;
    s(idx).match_status = 'Unmatched CSV';
    s(idx).timeBase = timeBaseSource;
    s(idx).startTime = NaN;
    s(idx).endTime = NaN;
    s(idx).duration_measured_s = NaN;
    s(idx).duration_param_s = csv_row.duration_ms / 1000;
    s(idx).frequency_measured = NaN;
    s(idx).frequency_param = csv_row.frequency_param;
    s(idx).stimType = csv_row.stimType{1};
    s(idx).contrast1 = csv_row.contrast1;
    s(idx).contrast2 = csv_row.contrast2;
    s(idx).time = [];
    s(idx).stimCHANGE = [];
end