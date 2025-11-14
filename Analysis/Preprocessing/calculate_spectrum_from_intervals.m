function [S, f, R, Serr] = calculate_spectrum_from_intervals(unit_spike_times, intervals, params)
% CALCULATE_SPECTRUM_FROM_INTERVALS computes the power spectrum from spike times.
%
%   Handles cases with empty intervals or no spikes by returning
%   S=[], f=[], R=0, Serr=[].

    % --- 1. Input Validation ---
    
    % Check for required 'params' input
    if nargin < 3
        error('Missing required input: params. Usage: [S,f,R,Serr] = calculate_spectrum_from_intervals(spikes, intervals, params)');
    end
    
    % Check that intervals, if not empty, is Nx2
    if ~isempty(intervals) && size(intervals, 2) ~= 2
        error('Input ''intervals'' must be an Nx2 matrix of [start, end] times.');
    end
    
    % --- 2. Data Preparation ---
    num_trials = size(intervals, 1);
    data = repmat(struct('spike_times', []), num_trials, 1); % init

    for i = 1:num_trials
        % Get the start and end time for the current trial
        startTime = intervals(i, 1);
        endTime   = intervals(i, 2);
        
        % Find logical indices of spikes within this trial (inclusive)
        spike_indices = (unit_spike_times >= startTime) & (unit_spike_times <= endTime);
        
        % Get spike times, subtract startTime, and convert to seconds
        trial_spikes = (unit_spike_times(spike_indices) - startTime) / 1000;
        
        % Store the zeroed spike times in our struct array
        data(i).spike_times = trial_spikes;
    end

    % --- 3. Handle "No Spikes" or "No Trials" Case ---
    
    % Consolidate all spikes found across all trials
    % 'vertcat' is a fast way to join all 'spike_times' fields
    all_spikes = vertcat(data.spike_times); 
    
    % --- 4. Calculate Spectrum  ---
    
    % If we get here, it means we have at least one spike
    [S, f, R, Serr] = mtspectrumpt(data, params, [], params.t);
    
end