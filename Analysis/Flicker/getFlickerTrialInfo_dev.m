Subject = 'M25065';
Session = '20250528';
AcquisitionsID = '0';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);

% get 1hz sync pulse files
nidq_1hz_file = dir(fullfile(nidqDir, '*nidq.xd_8_3_500.txt'));
nidq_tprime_file = dir(fullfile(nidqDir, 'nidq_sync_tprime.txt'));

% Build the full file paths
nidq_1hz_file = fullfile(nidq_1hz_file.folder, nidq_1hz_file.name);
nidq_tprime_file = fullfile(nidq_tprime_file.folder, nidq_tprime_file.name);

%% get nidq data for this gX recording

igkey=0; % % choose which stimulus to extract info from
nidq_channels.stimON = 5; % which channels to extract and their names
nidq_channels.stimCHANGE = 6;
nidq_channels.asyncPulse = 1;

doPlot=true;

% get minimally processed nidq data
nidq = extractNidqData(ExpInfo, igkey, nidqDir, nidq_channels, nidq_1hz_file, nidq_tprime_file,doPlot);

%% generate trials struct from nidq and trialParams .csv file

g_idx = find(ExpInfo.gKey == igkey);
this_trialParamsCSVfile = ExpInfo.TrialParamsFiles{g_idx};

trials = crateTrialsStruct_flicker(nidq, this_trialParamsCSVfile)

%%

% get sync pulse times
thisSync = NidqSyncChan(nidqidx);
threshold = mean(thisSync);
thisSync(thisSync<threshold)=0;
thisSync(thisSync>=threshold)=1;
nidq_sync_idx = find(diff(thisSync)==1);
syncTimes_nidq = thisTime(nidq_sync_idx+1);

%% wheel analysis

wheel_tbl = readtable(ExpInfo.WheelFiles{stimIdx});
wheel_tbl.dist = wheel2unit(wheel_tbl.Wheel,4096, 19.5);
wheel_tbl.ddist_dt = [nan; diff(wheel_tbl.dist)./diff(wheel_tbl.ArduinoTime/1000)];
wheel_tbl.ArduinoTime_mid = movmean(wheel_tbl.ArduinoTime/1000,2);
wheel_tbl(1,:)=[];

syncTimes_bonsai = unique(wheel_tbl.LastSyncPulseTime)/1000;

wheel_tbl.nidqTime = mapTimestampsUsingAsyncPulse(syncTimes_bonsai,syncTimes_nidq,wheel_tbl.ArduinoTime_mid);

wheel_tbl.time = wheel_tbl.nidqTime;
wheel_tbl.smthSpeed = wheel_tbl.ddist_dt;



%% get stimulus start times for flicker

stim_idx_to_use = stimIdx;

ExpInfo.Stimulus(stim_idx_to_use(1))
allStimTimes = [];
allStimTimesImec0=[];



    if stimIdx<height(ExpInfo)
        idx = [ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx+1)];
    else
        idx = [ExpInfo.smp_nidq(stimIdx):numel(NidqTime)];
    end



    % get PD signal and timestamps
    thisPD = stimON(idx);
    thisTime = NidqTime(idx);
    stimTimeIdx = idx;

    % convert to a digital signal
    PD_binary = thisPD;
    threshold = mean(stimON);
    PD_binary(PD_binary<threshold)=0;
    PD_binary(PD_binary>=threshold)=1;

    % find trial start times
    idx = find(diff(PD_binary)==1); % indexes of where PD goes up
    idx = idx+1;
    diffTimes = thisTime(idx); % time stamps of where PD goes up
    diffDiffTimes = diff(diffTimes); % time between PD going up.
    idx2 = find(diffDiffTimes>1.5); % For new trials its ~2, i.e. >1.5s
    idx2=idx2+1;
    trialStartIdx = idx(idx2); % convert back to original PD index
    trialStartIdx = [idx(1), trialStartIdx]; % add the first trial
    stimTimes = thisTime(trialStartIdx); % get the PD timestamp for trial start idx

    timeFrame = [NidqTime(stimTimeIdx(1)), NidqTime(stimTimeIdx(end))-5];

    % convert stim times to imec0 timestamp using interp1
    % (maybe this should be converted to tprime)
    this_nidqsync_nidq = nidq_syncTimes_nidq(nidq_syncTimes_nidq >= timeFrame(1) & nidq_syncTimes_nidq <= timeFrame(2));
    this_nidqsync_imec0 = nidq_syncTimes_imec0(nidq_syncTimes_imec0 >= timeFrame(1) & nidq_syncTimes_imec0 <= timeFrame(2));
    this_nidqsync_nidq=this_nidqsync_nidq(:); this_nidqsync_imec0=this_nidqsync_imec0(:);


    while numel(this_nidqsync_imec0)>numel(this_nidqsync_nidq)
        dv = this_nidqsync_imec0(1:numel(this_nidqsync_nidq))-this_nidqsync_nidq;
        idx = find(abs(dv)>0.5,1,'first');
        this_nidqsync_imec0(idx)=[];
    end

    while numel(this_nidqsync_nidq)>numel(this_nidqsync_imec0)
        dv = this_nidqsync_nidq(1:numel(this_nidqsync_imec0))-this_nidqsync_imec0;
        idx = find(abs(dv)>0.5,1,'first');
        this_nidqsync_nidq(idx)=[];
    end


    stimTimes_imec0 = interp1(this_nidqsync_nidq, this_nidqsync_imec0, stimTimes);

    allStimTimes = cat(2,allStimTimes, stimTimes);
    % allStimTimesImec0 = cat(2,allStimTimesImec0, stimTimes_imec0);





% get trial params
% quick hack for fake photodiode signal
photodiode.up_times = allStimTimes;
photodiode.down_times = photodiode.up_times+2;

allFlickerTrials = createTrialsStruct(ExpInfo,stim_idx_to_use,photodiode,wheel_tbl);


% rename struct fields
% Define the new field names in the desired order
new_fields = {'flicker_type', 'duration', 'frequency', 'phase_1', 'phase_2', 'contrast_1', 'contrast_2', 'bonsai_time', 'arduino_time', 'start_time', 'stop_time', 'WheelSpeedTrace', 'meanWheelSpeed' };
% Convert the struct to a cell array
data_as_cells = struct2cell(allFlickerTrials);
% Recreate the struct with the new field names
allFlickerTrials = cell2struct(data_as_cells, new_fields, 1);


stat_idx = [allFlickerTrials.meanWheelSpeed]<=1;
run_idx = [allFlickerTrials.meanWheelSpeed]>3;
[allFlickerTrials.runFlag] = deal([nan]);
[allFlickerTrials(stat_idx).runFlag] = deal([0]);
[allFlickerTrials(run_idx).runFlag] = deal([1]);


% get set of intervals for trials
statTrials = allFlickerTrials([allFlickerTrials.runFlag]==0);
runTrials = allFlickerTrials([allFlickerTrials.runFlag]==1);

flickFreqs = unique([allFlickerTrials.frequency]);



%% chronux analysis - green

%% chronux analysis - uv

params.Fs=1000;
params.fpass = [0.5 100];
params.trialave = 0;
params.pad = 0;
params.tapers = [2 4];
params.err =[1 0.05];

alllocs=cat(1,units.location);
idx = find(alllocs(:,2)<4000);

for ii = 1:numel(idx)

    iunit = idx(ii);

    try
figure
subplot(451)


    
% spikeTimes = cat(1,units.spike_times).*1000;
spikeTimes = units(iunit).spike_times.*1000;


for ifreq = 1:numel(flickFreqs)
    % stationary/locomotion trials
    stimTrials = statTrials([statTrials.frequency]==flickFreqs(ifreq) & [statTrials.contrast_1]==1);
    stat_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;

    stimTrials = runTrials([runTrials.frequency]==flickFreqs(ifreq) & [runTrials.contrast_1]==1);
    run_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;
end



     for ifreq = 1:numel(flickFreqs);

    tints = stat_intervalsCell{ifreq};
    trialLength = round(mean(tints(:,2)-tints(:,1)),-1,'decimals'); % round to 10ms

    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    intStarts = stimOnTimes;
    intStops = stimOnTimes+trialLength;

    ipsth = 1;
    % loop through intervals (usually trials) and bin spikes
    clear p stat_t

    for itrial = 1:size(tints,1)

        % spikes that occur during this interval
        relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
            spikeTimes <= intStops(itrial));
        % 0-centre to stimulus onset. save these for later raster plot
        p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
        p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
        stat_t(itrial).spike_times = p(1).trial(itrial).trialCentricSpikes/1000;

    end


    [S,f,R,Serr]=mtspectrumpt(stat_t,params);
    stat_t(itrial).S = S; stat_t(itrial).f=f; stat_t(itrial).R=R; stat_t(itrial).Serr = Serr;
    stat_S=S; stat_f = f; stat_R = R; stat_Serr = Serr;
    % stat_S = pow2db(stat_S);

    % run
    % 
    tints = run_intervalsCell{ifreq};
    trialLength = round(mean(tints(:,2)-tints(:,1)),-1,'decimals'); % round to 10ms

    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    intStarts = stimOnTimes;
    intStops = stimOnTimes+trialLength;

    ipsth = 1;
    % loop through intervals (usually trials) and bin spikes
    clear p

    for itrial = 1:size(tints,1)

        % spikes that occur during this interval
        relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
            spikeTimes <= intStops(itrial));
        % 0-centre to stimulus onset. save these for later raster plot
        p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
        p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
        run_t(itrial).spike_times = p(1).trial(itrial).trialCentricSpikes/1000;

    end

    [S,f,R,Serr]=mtspectrumpt(run_t,params);

    run_t(itrial).S = S; run_t(itrial).f=f; run_t(itrial).R=R; run_t(itrial).Serr = Serr;
    run_S=S; run_f = f; run_R = R; run_Serr = Serr;
    % run_S = pow2db(run_S);

    subplot(4,5,ifreq)
    shadedErrorBar(f,nanmean(stat_S,2),nansem(stat_S,2))
    shadedErrorBar(f,nanmean(run_S,2),nansem(run_S,2),'lineprops','r')
    ax=gca; ax.XTick = [10:10:100]; grid on
     % ylim([30 50])
    title(flickFreqs(ifreq))
end

% Get all subplot axes from the current figure (excluding legends)
allAxes = findall(gcf, 'type', 'axes');
isLegend = arrayfun(@(x) strcmp(get(x, 'Tag'), 'legend'), allAxes);
allAxes(isLegend) = [];

% Initialize global min and max with extreme values
global_ymin = inf;
global_ymax = -inf;

% First loop: Find the widest y-axis range among all subplots
for i = 1:length(allAxes)
    current_ylim = get(allAxes(i), 'YLim');
    global_ymin = min(global_ymin, current_ylim(1));
    global_ymax = max(global_ymax, current_ylim(2));
end

% Second loop: Apply the determined global limits to all subplots
for i = 1:length(allAxes)
    set(allAxes(i), 'YLim', [global_ymin, global_ymax]);
end
 % pause
 % close
sgtitle('GREEN, population spiking')


% UV

for ifreq = 1:numel(flickFreqs)
    % stationary/locomotion trials
    stimTrials = statTrials([statTrials.frequency]==flickFreqs(ifreq) & [statTrials.contrast_2]==1);
    stat_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;

    stimTrials = runTrials([runTrials.frequency]==flickFreqs(ifreq) & [runTrials.contrast_2]==1);
    run_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;
end

     for ifreq = 1:numel(flickFreqs);

    tints = stat_intervalsCell{ifreq};
    trialLength = round(mean(tints(:,2)-tints(:,1)),-1,'decimals'); % round to 10ms

    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    intStarts = stimOnTimes;
    intStops = stimOnTimes+trialLength;

    ipsth = 1;
    % loop through intervals (usually trials) and bin spikes
    clear p stat_t

    for itrial = 1:size(tints,1)

        % spikes that occur during this interval
        relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
            spikeTimes <= intStops(itrial));
        % 0-centre to stimulus onset. save these for later raster plot
        p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
        p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
        stat_t(itrial).spike_times = p(1).trial(itrial).trialCentricSpikes/1000;

    end


    [S,f,R,Serr]=mtspectrumpt(stat_t,params);
    stat_t(itrial).S = S; stat_t(itrial).f=f; stat_t(itrial).R=R; stat_t(itrial).Serr = Serr;
    stat_S=S; stat_f = f; stat_R = R; stat_Serr = Serr;
    % stat_S = pow2db(stat_S);

    % run
    % 
    tints = run_intervalsCell{ifreq};
    trialLength = round(mean(tints(:,2)-tints(:,1)),-1,'decimals'); % round to 10ms

    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    intStarts = stimOnTimes;
    intStops = stimOnTimes+trialLength;

    ipsth = 1;
    % loop through intervals (usually trials) and bin spikes
    clear p

    for itrial = 1:size(tints,1)

        % spikes that occur during this interval
        relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
            spikeTimes <= intStops(itrial));
        % 0-centre to stimulus onset. save these for later raster plot
        p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
        p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
        run_t(itrial).spike_times = p(1).trial(itrial).trialCentricSpikes/1000;

    end

    [S,f,R,Serr]=mtspectrumpt(run_t,params);

    run_t(itrial).S = S; run_t(itrial).f=f; run_t(itrial).R=R; run_t(itrial).Serr = Serr;
    run_S=S; run_f = f; run_R = R; run_Serr = Serr;
    % run_S = pow2db(run_S);

    subplot(4,5,ifreq+10)
    shadedErrorBar(f,nanmean(stat_S,2),nansem(stat_S,2))
    shadedErrorBar(f,nanmean(run_S,2),nansem(run_S,2),'lineprops','r')
    ax=gca; ax.XTick = [10:10:100]; grid on
     % ylim([30 50])
    title(flickFreqs(ifreq))
end

% Get all subplot axes from the current figure (excluding legends)
allAxes = findall(gcf, 'type', 'axes');
isLegend = arrayfun(@(x) strcmp(get(x, 'Tag'), 'legend'), allAxes);
allAxes(isLegend) = [];

% Initialize global min and max with extreme values
global_ymin = inf;
global_ymax = -inf;

% First loop: Find the widest y-axis range among all subplots
for i = 1:length(allAxes)
    current_ylim = get(allAxes(i), 'YLim');
    global_ymin = min(global_ymin, current_ylim(1));
    global_ymax = max(global_ymax, current_ylim(2));
end

% Second loop: Apply the determined global limits to all subplots
for i = 1:length(allAxes)
    set(allAxes(i), 'YLim', [global_ymin, global_ymax]);
end
 % pause
 % close

sgtitle(sprintf('Unit: %d', iunit))

pause

close all

    catch
        disp(sprintf('Skipping unit %d', iunit));
        close all 

    end

end


%%


options.intervalStart = 0;
options.intervalEnd = 2;
options.binSpacing=2;
options.uniqueVals = {[2 5 10 20 30 40 50 60 70 200], [0 1], [0 1]};

[anUnits, cond] = getBinnedSpikeCounts(allFlickerTrials, units, {'frequency','contrast_1','runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

for iunit = 1:numel(units)
    units(iunit).tuning_mean = cellfun(@mean, units(iunit).allSpikes,'UniformOutput',false);
    units(iunit).tuning_sem = cellfun(@(x) sem(x,2), units(iunit).allSpikes,'UniformOutput',false);
end


%% plot tuning curves

for iunit = 1:numel(units)
    subplot(211)
    title('green')
    shadedErrorBar(1:10, [units(iunit).tuning_mean{:,2,1}], [units(iunit).tuning_sem{:,2,1}], 'lineProps','k')
    shadedErrorBar(1:10, [units(iunit).tuning_mean{:,2,2}], [units(iunit).tuning_sem{:,2,2}], 'lineProps','r')
   ax=gca; ax.XTick=1:10; ax.XTickLabels = flickFreqs
    defaultAxesProperties(gca, true)
    subplot(212)
    title('uv')
    shadedErrorBar(1:10, [units(iunit).tuning_mean{:,1,1}], [units(iunit).tuning_sem{:,1,1}],'lineProps', 'k')
    shadedErrorBar(1:10, [units(iunit).tuning_mean{:,1,2}], [units(iunit).tuning_sem{:,1,2}],'lineProps', 'r')
   ax=gca; ax.XTick=1:10; ax.XTickLabels = flickFreqs
    defaultAxesProperties(gca, true)
    sgtitle(iunit)
    pause
close
end


%% plot mean tuning curves
% --- Main Code to Calculate Mean and Standard Error ---
% This script assumes your struct array is named 'units' and is in the workspace.

% Get the total number of units (structs) in the array
num_units = numel(units);

% Check if there are enough units to calculate statistics
if num_units < 2
    error('Need at least two units to calculate standard error.');
end

% Get the dimensions from the first unit's 'tuning_mean' field.
% These correspond to: (frequency, colour, state)
[freq_bins, colour_size, state_size] = size(units(1).tuning_mean);

% Step 1: Aggregate all data into a single 4D matrix.
% Dimensions will be: (frequency, colour, state, unit_index)
all_data_matrix = zeros(freq_bins, colour_size, state_size, num_units);

for i = 1:num_units
    if isfield(units(i), 'tuning_mean') && iscell(units(i).tuning_mean)
        % Convert the cell array for the current unit to a matrix and store it
        all_data_matrix(:, :, :, i) = cell2mat(units(i).tuning_mean);
    else
        warning('Unit %d does not contain a valid ''tuning_mean'' field. Filling with NaN.', i);
        all_data_matrix(:, :, :, i) = NaN;
    end
end

% Step 2: Calculate the mean across the 4th dimension (across all units)
mean_tuning_curves = mean(all_data_matrix, 4, 'omitnan');

% Step 3: Calculate the standard deviation across the 4th dimension
std_dev_tuning_curves = std(all_data_matrix, 0, 4, 'omitnan');

% Step 4: Calculate the standard error of the mean (SEM)
% This correctly handles any potential missing data (NaNs)
num_valid_units = sum(~isnan(all_data_matrix), 4);
sem_tuning_curves = std_dev_tuning_curves ./ sqrt(num_valid_units);


% The 'mean_tuning_curves' and 'sem_tuning_curves' are now 10x2x2 matrices.

% --- Accessing the Results by Condition Name ---

% Example: Get the mean and SEM for the first colour and first state
mean_colour1_state1 = mean_tuning_curves(:, 1, 1);
sem_colour1_state1  = sem_tuning_curves(:, 1, 1);

% Example: Get the mean and SEM for the first colour and second state
mean_colour1_state2 = mean_tuning_curves(:, 1, 2);
sem_colour1_state2  = sem_tuning_curves(:, 1, 2);

% Example: Get the mean and SEM for the second colour and first state
mean_colour2_state1 = mean_tuning_curves(:, 2, 1);
sem_colour2_state1  = sem_tuning_curves(:, 2, 1);

% Example: Get the mean and SEM for the second colour and second state
mean_colour2_state2 = mean_tuning_curves(:, 2, 2);
sem_colour2_state2  = sem_tuning_curves(:, 2, 2);


% --- Optional: Visualize the Results with Error Bars ---

figure;
set(gcf, 'color', 'w'); % Set figure background to white

% Create an x-axis for the frequency bins (e.g., 1, 2, 3...)
frequency_axis = 1:freq_bins;

% Plot for State 1
subplot(1, 2, 1); % Create a subplot for the first state
hold on;
errorbar(frequency_axis, mean_colour1_state1, sem_colour1_state1, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Colour 1');
errorbar(frequency_axis, mean_colour2_state1, sem_colour2_state1, 's-', 'LineWidth', 1.5, 'DisplayName', 'Colour 2');
hold off;
title('State 1');
xlabel('Frequency Bin');
ylabel('Mean Response');
legend('show', 'Location', 'best');
grid on;
xlim([0.5, freq_bins + 0.5]); % Adjust x-axis for better appearance

% Plot for State 2
subplot(1, 2, 2); % Create a subplot for the second state
hold on;
errorbar(frequency_axis, mean_colour1_state2, sem_colour1_state2, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Colour 1');
errorbar(frequency_axis, mean_colour2_state2, sem_colour2_state2, 's-', 'LineWidth', 1.5, 'DisplayName', 'Colour 2');
hold off;
title('State 2');
xlabel('Frequency Bin');
ylabel('Mean Response');
legend('show', 'Location', 'best');
grid on;
xlim([0.5, freq_bins + 0.5]);

sgtitle('Mean Frequency Tuning Curves by State and Colour'); % Add a main title to the figure


%%

subplot(121)
shadedErrorBar(frequency_axis, mean_colour2_state1, sem_colour2_state1, 'lineProps', 'k')
shadedErrorBar(frequency_axis, mean_colour2_state2, sem_colour2_state2, 'lineProps', 'r')
title('green')

subplot(122)
shadedErrorBar(frequency_axis, mean_colour1_state1, sem_colour1_state1, 'lineProps', 'k')
shadedErrorBar(frequency_axis, mean_colour1_state2, sem_colour1_state2, 'lineProps', 'r')
title('uv')
