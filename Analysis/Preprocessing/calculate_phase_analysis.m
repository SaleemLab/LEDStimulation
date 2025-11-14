function [unit_phase_analysis] = calculate_phase_analysis(units_struct, trials_struct, config)
% CALCULATE_PHASE_ANALYSIS Performs batch phase-locking analysis for all units.
%
%   [unit_phase_analysis] = calculate_phase_analysis(units_struct, trials_struct, config)
%
%   This function loops through all trials ONCE and assigns spikes from all
%   units, which is much more efficient than looping through units first.
%
%   INPUTS:
%   - units_struct:  The main 'units' struct array (e.g., units(i).spike_times)
%   - trials_struct: The main 'trials' struct array (e.g., trials(j).startTime)
%   - config:        A struct with analysis parameters and field names.
%       config.uniqueFreqs    = [2, 5, 10, ...];
%       config.uniqueLums     = [1, 2];
%       config.uniqueColors   = [-1, 1];
%       config.uniqueStates   = [0, 1];
%       config.field_luminance = 'luminance';  % Field name in trials_struct
%       config.field_color     = 'contrast_coded';
%       config.field_frequency = 'frequency';
%       config.field_state     = 'is_running';
%
%   OUTPUT:
%   - unit_phase_analysis: A new struct array, 1-per-unit, containing
%     all spike phases and calculated vector strengths for all conditions.
%

% --- 1. Initialization ---

% Get condition dimensions from config
nUnits = numel(units_struct);
nFreqs = numel(config.uniqueFreqs);
nLums = numel(config.uniqueLums);
nColors = numel(config.uniqueColors);
nStates = numel(config.uniqueStates);
nTrials = numel(trials_struct);

% Define a struct for our conditions
condition_template = struct(...
    'all_spike_phases', [], ...    % Pooled spikes from all trials
    'trial_spike_phases', {{}}, ... % Cell array, 1 cell per trial
    'vector_strength', NaN, ...
    'mean_phase_rad', NaN, ...
    'p_value', NaN, ...
    'n_spikes', 0);

% Pre-allocate the entire output struct array
unit_phase_analysis(nUnits) = struct();
for iunit = 1:nUnits
    unit_phase_analysis(iunit).unit_id = units_struct(iunit).cluster_id; % Store unit ID
    unit_phase_analysis(iunit).conditions = repmat(condition_template, ...
        nColors, nLums, nFreqs, nStates); 
end

fprintf('Pre-allocation complete. Looping through %d trials...\n', nTrials);

% --- 2. Main Loop: Collect Spike Phases (Trial-centric) ---

for iTrial = 1:nTrials
    
    if mod(iTrial, 500) == 0
        fprintf('  Processing trial %d of %d\n', iTrial, nTrials);
    end
    
    % --- A. Get trial info ---
    current_trial = trials_struct(iTrial);
    t_start = current_trial.startTime;
    t_end = current_trial.endTime;

    % --- B. Get condition indices ---
    try
        c_idx = find(config.uniqueColors == current_trial.(config.field_color));
        l_idx = find(config.uniqueLums == current_trial.(config.field_luminance));
        f_idx = find(config.uniqueFreqs == current_trial.(config.field_frequency));
        s_idx = find(config.uniqueStates == current_trial.(config.field_state));
    catch e
        fprintf('Warning: Could not find field name in trials struct (e.g., %s). Skipping trial %d.\n', e.message, iTrial);
        continue;
    end
    
    % Skip if this trial is a control (e.g., 200Hz) or condition not found
    if isempty(c_idx) || isempty(l_idx) || isempty(f_idx) || isempty(s_idx)
        continue;
    end
    
    % --- C. Get correct phase lookup table ---
    trial_time_vec = current_trial.time_1khz;
    if isempty(trial_time_vec)
        continue; % Skip trial if time vector is missing
    end
    
    if current_trial.(config.field_color) == 1 % Green (adjust logic if needed)
        trial_phase_vec = current_trial.phase_green_rad;
    elseif current_trial.(config.field_color) == -1 % UV
        trial_phase_vec = current_trial.phase_uv_rad;
    else
        continue; % Skip if color code is not recognized
    end

    % --- D. Loop through all units for this trial ---
    for iunit = 1:nUnits
        
        % Find all spikes from this unit that occurred in this trial
        unit_spikes = units_struct(iunit).spike_times;
        trial_spike_times = unit_spikes(unit_spikes >= t_start & unit_spikes <= t_end);
        
        if isempty(trial_spike_times)
            continue; % No spikes from this unit in this trial
        end
        
        % --- E. Look up phase for each spike ---
        trial_spike_phases = zeros(size(trial_spike_times)); % Pre-allocate
        
        for iSpike = 1:numel(trial_spike_times)
            spike_time = trial_spike_times(iSpike);
            
            % Find the index of the closest time in the 1kHz time vector
            [~, time_idx] = min(abs(trial_time_vec - spike_time));
            
            % Get the phase at that exact time
            trial_spike_phases(iSpike) = trial_phase_vec(time_idx);
        end
        
        % --- F. Store the spike phases ---
        
        % Append this trial's spike phases to the correct condition bin
        unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).all_spike_phases = ...
            [unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).all_spike_phases; ...
             trial_spike_phases];

        % Append the vector of phases as a new cell in the trial_spike_phases field
        unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).trial_spike_phases{end+1} = ...
             trial_spike_phases;
    end
end

fprintf('All spike phases collected. Now calculating vector strengths...\n');

% --- 3. Calculate Vector Strength for all conditions ---

for iunit = 1:nUnits
    for c_idx = 1:nColors
        for l_idx = 1:nLums
            for f_idx = 1:nFreqs
                for s_idx = 1:nStates
                    
                    % Get all the spikes we collected for this condition
                    all_phases = unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).all_spike_phases;
                    n_spikes = numel(all_phases);
                    
                    unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).n_spikes = n_spikes;
                    
                    if n_spikes > 5 % Need a few spikes to calculate
                        
                        % --- This is the core calculation ---
                        % Convert phases (0 to 2pi) to complex vectors on the unit circle
                        complex_vectors = exp(1i * all_phases);
                        
                        % The mean vector is the average of all individual spike vectors
                        mean_complex_vector = mean(complex_vectors);
                        
                        % The Vector Strength is the length (magnitude) of the mean vector
                        vs = abs(mean_complex_vector);
                        
                        % The Mean Phase is the angle of the mean vector
                        phase = angle(mean_complex_vector);
                        
                        % Store the results
                        unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).vector_strength = vs;
                        unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).mean_phase_rad = phase;
                        
                        % Calculate significance (Requires Circular Stats Toolbox)
                        try
                            p_val = circ_rtest(all_phases);
                            unit_phase_analysis(iunit).conditions(c_idx, l_idx, f_idx, s_idx).p_value = p_val;
                        catch
                            % circ_rtest not found, p_value remains NaN
                        end
                    end
                end
            end
        end
    end
end

fprintf('Batch phase analysis complete.\n');

end

