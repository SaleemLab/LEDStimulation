%% Sinusoidal Flicker - single unit analysis

% 1. PSTH Plots: Show the trial-averaged PSTH to visualize the shape and timing of the response.
% 2. Mean Firing Rate (F0) Tuning Curve: Plot the mean firing rate during the stimulus. Show the pre-stimulus baseline rate as a horizontal line for comparison.
% 3. Normalized Power (F1) Tuning Curve: Plot the normalized power at the stimulus frequency, P_normalized(F1). This shows the strength of stimulus entrainment.
% 4. Vector Strength Tuning Curve: The primary measure of temporal fidelity.
% 5. Mean Phase Tuning Curve: Shows the response latency. Plot the mean phase angle of the response, which often reveals a systematic phase lag with increasing frequency.
% 6. Non-Linearity Tuning Curves: Plot the ratios for the first harmonic (P_norm(2F1)/P_norm(F1)) and first subharmonic (P_norm(F1/2)/P_norm(F1)).

%% load preprocessed session data

load('M25065_20250528_flickerSession.mat')

%% split units according to location

unit_depth_thresh=4000; % thalamus vs cortex
all_locs = cat(1,units.location);
lgn_indices = all_locs(:,2) < unit_depth_thresh;
cortex_indices = all_locs(:,2) >= unit_depth_thresh;

% Assign the 'area' field directly
[units(lgn_indices).area] = deal({'thalamus'});
[units(cortex_indices).area] = deal({'cortex'});

%% Basic processing of trials

% assign luminance parameter to trials
[trials(arrayfun(@(t) ismember('10^3', t.properties), trials)).luminance] = deal(3);
[trials(arrayfun(@(t) ismember('10^2', t.properties), trials)).luminance] = deal(2);

% add coded contrast field (useful when contrast_green/uv are redundant)
coded_cells = num2cell([trials.contrast_green] - [trials.contrast_uv]);
[trials.contrast_coded] = coded_cells{:};

% assign behavioural state to each trial
run_idx = find(cellfun(@(x) prop(x>0.5)>=0 & mean(x)>1, {trials.wheelSpeed_stim_period}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0 & mean(x)<1, {trials.wheelSpeed_stim_period}));

[trials.runFlag] = deal([nan]);
[trials(stat_idx).runFlag] = deal([0]);
[trials(run_idx).runFlag] = deal([1]);


%% get start/stop intervals of trials for PSTHs and power spectrums

% Define the conditions to loop through
luminance_vals = [2 3];
contrast_vals = [-1 1]; % -1 for green, 1 for uv
contrast_labels = {'Green', 'UV'};

% Create a cell array to hold the interval data for each condition
% Dimensions: {contrast, luminance}
allIntervalData = cell(numel(contrast_vals), numel(luminance_vals));

% Loop over all conditions to prepare the data
for iContrast = 1:numel(contrast_vals)
    for iLum = 1:numel(luminance_vals)

        conditionData = struct(); % Struct to hold data for this specific condition

        % Select trials for the current condition
        current_trials = trials([trials.contrast_coded] == contrast_vals(iContrast) & ...
            [trials.luminance] == luminance_vals(iLum));

        if isempty(current_trials)
            allIntervalData{iContrast, iLum} = []; % Mark as empty if no trials
            continue; % Skip to next condition
        end

        % Separate into stationary and running trials
        stat_trials = current_trials([current_trials.runFlag] == 0);
        run_trials  = current_trials([current_trials.runFlag] == 1);
        uniqueFreqs = unique([current_trials.frequency_param]);

        % Prepare intervals for stimuli (all freqs except baseline)
        stat_intervalsCell = {};
        run_intervalsCell = {};
        for ifreq = 1:numel(uniqueFreqs)-1
            % Stationary
            stimTrials_stat = stat_trials([stat_trials.frequency_param] == uniqueFreqs(ifreq));
            if ~isempty(stimTrials_stat)
                stat_intervalsCell{ifreq} = [vertcat(stimTrials_stat.startTime), vertcat(stimTrials_stat.endTime)] .* 1000;
            else
                stat_intervalsCell{ifreq} = [];
            end
            % Running
            stimTrials_run = run_trials([run_trials.frequency_param] == uniqueFreqs(ifreq));
            if ~isempty(stimTrials_run)
                run_intervalsCell{ifreq} = [vertcat(stimTrials_run.startTime), vertcat(stimTrials_run.endTime)] .* 1000;
            else
                run_intervalsCell{ifreq} = [];
            end
        end

        % Prepare intervals for baseline (last frequency)
        ifreq_baseline = numel(uniqueFreqs);
        % Stationary
        stimTrials_stat_blank = stat_trials([stat_trials.frequency_param] == uniqueFreqs(ifreq_baseline));
        if ~isempty(stimTrials_stat_blank)
            stat_blankIntervals = [vertcat(stimTrials_stat_blank.startTime), vertcat(stimTrials_stat_blank.endTime)] .* 1000;
        else
            stat_blankIntervals = [];
        end
        % Running
        stimTrials_run_blank = run_trials([run_trials.frequency_param] == uniqueFreqs(ifreq_baseline));
        if ~isempty(stimTrials_run_blank)
            run_blankIntervals = [vertcat(stimTrials_run_blank.startTime), vertcat(stimTrials_run_blank.endTime)] .* 1000;
        else
            run_blankIntervals = [];
        end

        % Store the prepared data
        conditionData.stat_intervals = stat_intervalsCell;
        conditionData.run_intervals = run_intervalsCell;
        conditionData.stat_baseline = stat_blankIntervals;
        conditionData.run_baseline = run_blankIntervals;

        allIntervalData{iContrast, iLum} = conditionData;
        allIntervalData{iContrast, iLum}.contrast = contrast_vals(iContrast);
        allIntervalData{iContrast, iLum}.luminance = luminance_vals(iLum);
    end
end


%% PSTHs
% Options for the makePSTH function
options.binWidth = 10;
options.preTime = 1000;
options.postTime = 1000;
options.smoothType = 'gaussian';
options.smoothWidth = 175;
options.plot = true; % makePSTH will handle the plotting on the current axes
options.cols = inferno(9);
options.psthSize = 1;
options.getReliability = false;


plotFlag=true;
if plotFlag
    for iunit = 115%:numel(filtUnits)
        fprintf('Processing unit %d...\n', iunit);

        figure('Position', [50, 50, 800, 1000]);
        plot_row = 1;

        for iContrast = 1:numel(contrast_vals)
            contrast_label = contrast_labels{iContrast};

            for iLum = 1:numel(luminance_vals)
                lum_val = luminance_vals(iLum);

                % Retrieve the pre-calculated data for this condition
                conditionData = allIntervalData{iContrast, iLum};

                % Skip if this condition had no trials, but advance plot row
                if isempty(conditionData)
                    plot_row = plot_row + 1;
                    continue;
                end

                % Plot stationary PSTH (left column)
                subplot(4, 2, (plot_row-1)*2 + 1);
                makePSTH(units(iunit).spike_times * 1000, conditionData.stat_intervals, conditionData.stat_baseline, options);
                title(sprintf('Stationary: %s, Lum %d', contrast_label, lum_val));

                % Plot running PSTH (right column)
                subplot(4, 2, (plot_row-1)*2 + 2);
                makePSTH(units(iunit).spike_times * 1000, conditionData.run_intervals, conditionData.run_baseline, options);
                title(sprintf('Running: %s, Lum %d', contrast_label, lum_val));

                plot_row = plot_row + 1;
            end
        end

        sgtitle(sprintf('Unit %d PSTHs', iunit));
        % pause;
        % close all;
    end
end
%% Mean Firing Rate (F0) Temporal Frequency Tuning Curves

% conditions: luminance, chromatic contrast (UV or green), behavioural state (e.g. stat vs run)

options.intervalStart = 0;
options.intervalEnd = 2; % trial duration
options.binSpacing = 2; % just one bin for tuning curves
options.trialStartField = 'startTime';
options.spikeTimesField = 'spike_times';
options.allSpikesField = 'spikes_tf';

uniqueFreqs = unique([trials.frequency_param]);
paramsOfInterest = {'frequency_param', 'contrast_coded', 'luminance', 'runFlag'};
options.uniqueVals = {uniqueFreqs, [-1 1], [2, 3], [0 1]};

[units, cond, varInfo] = getBinnedSpikeCounts(trials, units, paramsOfInterest, options);

%% plot tuning curves across conditions

plotOptions = struct();

plotOptions.xAxisDimName = 'frequency_param';
plotOptions.subplotDimNames = {'contrast_coded', 'luminance'};
plotOptions.overlayDimNames = {'runFlag'};
plotOptions.baseline = 'last';
plotOptions.baselineCollapseDims = {'contrast_coded'};
plotOptions.plotStyles = {'k-', 'r-'};
contrastMap = containers.Map([-1, 1], {'Green', 'Ultraviolet'});
plotOptions.valueLabels.contrast_coded = contrastMap;
luminanceMap = containers.Map([2, 3], {'10^2', '10^3'});
plotOptions.valueLabels.luminance = luminanceMap;
stateMap = containers.Map([0, 1], {'Stat', 'Run'});
plotOptions.valueLabels.runFlag = stateMap;

plotOptions.plotIndividualPoints = false;

plotFlag=true;
if plotFlag
    % for iunit = 1:numel(units)

    iunit = 115
    plotFlickerTuningCurves(units(iunit).spikes_tf, varInfo, plotOptions)
    % pause
    % close
    % end
end

%% Spiking power spectrums
% power at F1, 2*F1 and F1/2, normalised by baseline (ratiometric)

% chronux mtspectrumpt options
params=struct;
params.Fs=30000;
params.fpass = [0.5 200];
params.trialave = 0;
params.pad = 1; % increase for finer f resolution
params.tapers = [2 4];
params.err =[1 0.05];

target_duration_ms = 2000; %trial duration
duration_tolerance_ms = 100;
tic
[units_power_analysis, f] = calculate_power_analysis(units, allIntervalData, uniqueFreqs, params, target_duration_ms, duration_tolerance_ms);
toc
%%
% get baseline data for each luminance by combining across colors

for iLum = 1:2
    stat_baseline_alllumi{iLum} = cat(1,allIntervalData{1,iLum}.stat_baseline, allIntervalData{2,iLum}.stat_baseline);
    run_baseline_alllumi{iLum} = cat(1,allIntervalData{1,iLum}.run_baseline, allIntervalData{2,iLum}.run_baseline);
end

% power spectrums for each stimulus frequency, for each condition

for iunit = 1:numel(units)
    unit_spike_times = units(iunit).spike_times .* 1000; % trial intervals are set to ms

    for iLum = 1:2

        % --- Stat Baseline ---
        intervals_stat = stat_baseline_alllumi{iLum};
        [S_stat, f_stat, R_stat, Serr_stat] = calculate_spectrum_from_intervals(unit_spike_times, intervals_stat, params);

        units(iunit).powerAnalysis.baselinePower_stat{iLum} = S_stat;
        units(iunit).powerAnalysis.freq_stat = f_stat;

        % --- Run Baseline ---
        intervals_run = run_baseline_alllumi{iLum};
        [S_run, f_run, R_run, Serr_run] = calculate_spectrum_from_intervals(unit_spike_times, intervals_run, params);

        units(iunit).powerAnalysis.baselinePower_run{iLum} = S_run;
        units(iunit).powerAnalysis.freq_run  = f_run;

        for iContrast = 1:2
            for ifreq = 1:numel(uniqueFreqs)-1

                % --- Stat Power Spectrum --- %
                intervals_stat = allIntervalData{iContrast,iLum}.stat_intervals{ifreq};
                [S_stat, f_stat, R_stat, Serr_stat] = calculate_spectrum_from_intervals(unit_spike_times, intervals_stat, params);

                units(iunit).powerAnalysis.power_stat{iContrast,iLum,ifreq} = S_stat;
                units(iunit).powerAnalysis.power_baselineNorm_stat{iContrast,iLum,ifreq} = ...
                    units(iunit).powerAnalysis.power_stat{iContrast,iLum,ifreq}./mean(units(iunit).powerAnalysis.baselinePower_stat{iLum},2);

                % power at F1, F1/2 and 2*F1
                f=f_stat;
                target = uniqueFreqs(ifreq); %F1
                [~, idx] = min(abs(f - target));                
                units(iunit).powerAnalysis.powerF1_stat{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_stat{iContrast,iLum,ifreq}(idx,:);
                units(iunit).powerAnalysis.powerF1norm_stat{iContrast,iLum,ifreq} = units(iunit).power_baselineNorm_stat{iContrast,iLum,ifreq}(idx,:);

                target  = uniqueFreqs(ifreq)*2; %F1*2
                [~, idx] = min(abs(f - target));
                units(iunit).powerAnalysis.powerF1_double_stat{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_stat{iContrast,iLum,ifreq}(idx,:);
                units(iunit).powerAnalysis.powerF1norm_double_stat{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_baselineNorm_stat{iContrast,iLum,ifreq}(idx,:);

                target  = uniqueFreqs(ifreq)/2; %F1/2
                [~, idx] = min(abs(f - target));
                units(iunit).powerAnalysis.powerF1_half_stat{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_stat{iContrast,iLum,ifreq}(idx,:);
                units(iunit).powerAnalysis.powerF1norm_half_stat{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_baselineNorm_stat{iContrast,iLum,ifreq}(idx,:);


                % --- Run Power Spectrum --- %
                intervals_run = allIntervalData{iContrast,iLum}.run_intervals{ifreq};
                [S_run, f_run, R_run, Serr_run] = calculate_spectrum_from_intervals(unit_spike_times, intervals_run, params);

                units(iunit).powerAnalysis.power_run{iContrast,iLum,ifreq} = S_run;
                units(iunit).powerAnalysis.power_baselineNorm_run{iContrast,iLum,ifreq} = ...
                    units(iunit).powerAnalysis.power_run{iContrast,iLum,ifreq}./mean(units(iunit).powerAnalysis.baselinePower_run{iLum},2);

                % power at F1, F1/2 and 2*F1
                f=f_run;
                target = uniqueFreqs(ifreq); %F1
                [~, idx] = min(abs(f - target));
                units(iunit).powerAnalysis.powerF1_run{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_run{iContrast,iLum,ifreq}(idx,:);
                units(iunit).powerAnalysis.powerF1norm_run{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_baselineNorm_run{iContrast,iLum,ifreq}(idx,:);

                target  = uniqueFreqs(ifreq)*2; %F1*2
                [~, idx] = min(abs(f - target));
                units(iunit).powerAnalysis.powerF1_double_run{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_run{iContrast,iLum,ifreq}(idx,:);
                units(iunit).powerAnalysis.powerF1norm_double_run{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_baselineNorm_run{iContrast,iLum,ifreq}(idx,:);

                target  = uniqueFreqs(ifreq)/2; %F1/2
                [~, idx] = min(abs(f - target));
                units(iunit).powerAnalysis.powerF1_half_run{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_run{iContrast,iLum,ifreq}(idx,:);
                units(iunit).powerAnalysis.powerF1norm_half_run{iContrast,iLum,ifreq} = units(iunit).powerAnalysis.power_baselineNorm_run{iContrast,iLum,ifreq}(idx,:);

            end

        end
    end



end


%% plot power spectrums for all frequencies (one stim condition) 
iContrast = 1; iLum = 2;

for iunit = 100:numel(units)
figure
for ifreq = 1:9
    ax(ifreq) = subplot(2,5,ifreq);
    S = units_power_analysis(iunit).powerAnalysis.power_stat{iContrast,iLum,ifreq};
    shadedErrorBar(f,nanmean(S,2),nansem(S,2))

    S = units_power_analysis(iunit).powerAnalysis.power_run{iContrast,iLum,ifreq};
    shadedErrorBar(f,nanmean(S,2),nansem(S,2), 'lineProps', 'r')
    % ax(ifreq).YTickLabel = log2(ax(ifreq).YTick);
    % refline(0,1)
    xlim([0 100])
    title(uniqueFreqs(ifreq))
    ax=gca; ax.XTick = [2 5 10:10:100]; grid on
    
    
end

linkaxes(ax,'y')

sgtitle(['Unit:', num2str(iunit), ', iContrast:', num2str(iContrast), ' iLum:', num2str(iLum)])

ifreq=10;
subplot(2,5,10);
S = units_power_analysis(iunit).powerAnalysis.baselinePower_stat{iLum};
shadedErrorBar(f,nanmean(S,2),nansem(S,2))
S = units_power_analysis(iunit).powerAnalysis.baselinePower_run{iLum};
shadedErrorBar(f,nanmean(S,2),nansem(S,2), 'lineProps', 'r')
title(200)
    xlim([0 100])
    ax=gca; ax.XTick = [2 5 10:10:100]; grid on

pause
close
end



%% plot F1 tuning curve

% Extract the cell array for the condition of interest
iContrast=1;
iLum=2;

data_cell_stat = units(115).powerAnalysis.powerF1norm_stat(iContrast,iLum,:);
data_cell_run = units(115).powerAnalysis.powerF1norm_run(iContrast,iLum,:);

figure;  hold on
errorbar(1:9, squeeze(cellfun(@(x) mean(x,2), data_cell_stat)), squeeze(cellfun(@(x) sem(x,2), data_cell_stat)), 'ko-', 'LineWidth', 1.5, 'MarkerSize', 8);
errorbar(1:9, squeeze(cellfun(@(x) mean(x,2), data_cell_run)), squeeze(cellfun(@(x) sem(x,2), data_cell_run)), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);

% Add labels for clarity
xlabel('Frequency Index');
ylabel('Normalized Power at F1 (Mean +/- SEM)');
title('Unit 115: F1 Power Tuning Curve');
grid on;

%% Plot F2/F1 tuning curve

% Extract the cell array for the condition of interest
iContrast=1;
iLum=2;

data_cellF1_stat = units(115).powerAnalysis.powerF1norm_stat(iContrast,iLum,:);
data_cellF1_run = units(115).powerAnalysis.powerF1norm_run(iContrast,iLum,:);

data_cellF2_stat = units(115).powerAnalysis.powerF1norm_double_stat(iContrast,iLum,:);
data_cellF2_run = units(115).powerAnalysis.powerF1norm_double_run(iContrast,iLum,:);

% Define a small threshold to avoid division by zero
threshold = 0.1; 

% Create a function handle that filters both vectors before dividing
% 1. Find indices 'idx' where f1 power is above the threshold
% 2. Use 'idx' to select elements from both f2 and f1
% 3. Perform the division
filter_and_divide = @(f2, f1) f2(f1 > threshold) ./ f1(f1 > threshold);

ratio_cell = cellfun(filter_and_divide, data_cellF2_stat, data_cellF1_stat, 'UniformOutput', false);
mean_ratios = cellfun(@(x) mean(x,2), ratio_cell);
sem_ratios  = cellfun(@(x) sem(x,2), ratio_cell); % Assumes you have a 'sem' function

% 2. Squeeze the 1x1x9 results into 9x1 vectors for plotting
y_means  = squeeze(mean_ratios);
y_errors = squeeze(sem_ratios);

% 3. Plot
figure;
errorbar(1:9, y_means, y_errors, 'ko-');
xlabel('Frequency Index');
ylabel('Harmonic Ratio (P_norm(2F1)/P_norm(F1))');
grid on;

hold on
ratio_cell = cellfun(filter_and_divide, data_cellF2_run, data_cellF1_run, 'UniformOutput', false);
mean_ratios = cellfun(@(x) mean(x,2), ratio_cell);
sem_ratios  = cellfun(@(x) sem(x,2), ratio_cell); % Assumes you have a 'sem' function
y_means  = squeeze(mean_ratios);
y_errors = squeeze(sem_ratios);
errorbar(1:9, y_means, y_errors, 'ro-');

%% phase analysis
config.uniqueFreqs    = uniqueFreqs;
config.uniqueLums     = [2 3];
config.uniqueColors   = [-1, 1]; % only full contrast tf trials
config.uniqueStates   = [0, 1]; % stat/run
config.field_luminance = 'luminance';  % Field name in trials_struct
config.field_color     = 'contrast_coded';
config.field_frequency = 'frequency_param';
config.field_state     = 'runFlag';

[unit_phase_analysis] = calculate_phase_analysis(units, trials, config)

%% phase analysis



%% polarplots of phase

for iunit = 372
% --- 1. Setup Constants ---
nbins = 45;
binedges = linspace(0, 2*pi, nbins + 1);
binvals = movmean(binedges, 2, 'EndPoints', 'discard'); % Get bin centers
binvals_closed = [binvals, binvals(1)]; % Add first bin center to the end

nFreqsToPlot = 4; % How many frequencies to overlay (e.g., first 3)
colors = inferno(nFreqsToPlot); % Colormap for the frequencies

% --- 2. Create the Figure ---
figure('Name', 'Stationary vs Running Phase Distributions', ...
       'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.5]);
sgtitle('Overlaid Spike Phase Distributions (Line Plot)', 'FontSize', 14);


% --- 3. Create Subplot for iState = 1 (Stationary) ---

% --- THIS IS THE FIX ---
% 1. Create a dummy subplot to get its position
temp_ax = subplot(1, 2, 1);
ax1_pos = temp_ax.Position;
% 2. Delete the dummy
delete(temp_ax);
% 3. Create the real polaraxes in that position
ax1 = polaraxes('Position', ax1_pos);
% --- END FIX ---

hold(ax1, 'on'); 
title(ax1, 'Stationary (iState = 1)');

iState = 1; % Index for Stationary
for ifreq = 1:nFreqsToPlot
    % Get the spike phases for this frequency and state
    % Note: Using (1,2,ifreq,iState) for (Color 1, Lum 2) as an example
    all_phases = unit_phase_analysis(iunit).conditions(1, 2, ifreq, iState).all_spike_phases;
    
    if ~isempty(all_phases)
        binCounts = histcounts(all_phases, binedges);
        binCounts_norm = binCounts / sum(binCounts); % Normalize to probability
        binCounts_closed = [binCounts_norm, binCounts_norm(1)]; % Close the loop
        
        polarplot(ax1, binvals_closed, binCounts_closed, ...
                  'LineWidth', 2, ...
                  'Color', [colors(ifreq, :), 0.7], ... % Add transparency
                  'DisplayName', sprintf('Freq %d', ifreq));
    end
end
legend(ax1, 'show');
hold(ax1, 'off');


% --- 4. Create Subplot for iState = 2 (Running) ---

% --- THIS IS THE FIX ---
% 1. Create a dummy subplot to get its position
temp_ax = subplot(1, 2, 2);
ax2_pos = temp_ax.Position;
% 2. Delete the dummy
delete(temp_ax);
% 3. Create the real polaraxes in that position
ax2 = polaraxes('Position', ax2_pos);
% --- END FIX ---

hold(ax2, 'on');
title(ax2, 'Running (iState = 2)');

iState = 2; % Index for Running
for ifreq = 1:nFreqsToPlot
    % Get the spike phases for this frequency and state
    all_phases = unit_phase_analysis(iunit).conditions(1, 2, ifreq, iState).all_spike_phases;
    
    if ~isempty(all_phases)
        binCounts = histcounts(all_phases, binedges);
        binCounts_norm = binCounts / sum(binCounts); % Normalize
        binCounts_closed = [binCounts_norm, binCounts_norm(1)]; % Close loop
        
        polarplot(ax2, binvals_closed, binCounts_closed, ...
                  'LineWidth', 2, ...
                  'Color', [colors(ifreq, :), 0.7], ...
                  'DisplayName', sprintf('Freq %d', ifreq));
    end
end
legend(ax2, 'show');
hold(ax2, 'off');


end