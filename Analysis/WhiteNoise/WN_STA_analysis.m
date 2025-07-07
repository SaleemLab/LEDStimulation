%% get session info

Subject = 'M25066';
Session = '20250701';
AcquisitionsID = '0';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);


%% load units

SorterName = 'kilosort4';
units = createClusterTable(ephysDir, SorterName);
units = table2struct(units);

%%
stimIdx=1; % choose the stimulus based on ExpInfo (gkey+1)
output = processStimulusData(stimIdx, ExpInfo, nidqDir);

%% load and process wheel data

wheel_tbl = readtable(ExpInfo.WheelFiles{stimIdx});
wheel_tbl.dist = wheel2unit(wheel_tbl.Wheel,4096, 19.5);
wheel_tbl.ddist_dt = [nan; diff(wheel_tbl.dist)./diff(wheel_tbl.ArduinoTime/1000)];
wheel_tbl.ArduinoTime_mid = movmean(wheel_tbl.ArduinoTime/1000,2);
wheel_tbl(1,:)=[];

syncTimes_bonsai = unique(wheel_tbl.LastSyncPulseTime)/1000;

wheel_tbl.nidqTime = mapTimestampsUsingAsyncPulse(syncTimes_bonsai,output.syncTimes_nidq,wheel_tbl.ArduinoTime_mid);

%% changepoints analysis to get stationary and locomotion epochs

wheelSpeed = wheel_tbl.ddist_dt;
wheel_zscoreSpeed = zscore(wheelSpeed);
wheelTime = wheel_tbl.nidqTime;
zThresh = 0.005; % moving standard deviations exceeded/fell below an empirical threshold of 0.005
sampleRate = 120; % of wheel data
smoothWin = 2; % moving standard deviation of speed (2s in Lohani)
timeBetween=0.5; % minimum time between off and the next on in seconds
minDuration=5; % min duration of state epoch
minLocoSpeed = 3;
maxStatSpeed = 1;
preBuffer =0.5; % buffer period to remove from start of epochs
postBuffer=0.5;

[stationary_intervals, locomotion_intervals] = ...
    changepointsWheelAnalysis(wheel_zscoreSpeed, wheelSpeed, wheelTime, zThresh, sampleRate,...
    smoothWin, timeBetween, minDuration, minLocoSpeed, maxStatSpeed,...
    preBuffer, postBuffer);


%% clear unnecessary variables

clear wheel_zscoreSpeed wheelSpeed wheelTime 


%% generate STAs

tic
nidqMetaFile = dir(fullfile(nidqDir, '*.meta'));
NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
niSampRate = str2double(NidqMeta.niSampRate);
timeWindowBefore = 0.2;
bufferVal = 0.5; % extra time buffer
nSampsBefore = ceil(timeWindowBefore * niSampRate);

for iunit = 1:numel(units)
    disp(['Processing Unit ', num2str(iunit)]);

    % 1. Extract and filter all spike times for the current unit
    spikeTimes = units(iunit).spike_times;
    spikeTimes = spikeTimes(spikeTimes > output.dutyCycleTime(1) + timeWindowBefore + bufferVal);
    spikeTimes = spikeTimes(spikeTimes < output.dutyCycleTime(end) - timeWindowBefore - bufferVal);
    
    % 2. Sort spikes into stationary and locomotion epochs
    is_stat = any((spikeTimes(:) >= stationary_intervals(:,1)') & (spikeTimes(:) <= stationary_intervals(:,2)'), 2);
    stat_SpikeTimes = spikeTimes(is_stat);

    is_loco = any((spikeTimes(:) >= locomotion_intervals(:,1)') & (spikeTimes(:) <= locomotion_intervals(:,2)'), 2);
    loco_SpikeTimes = spikeTimes(is_loco);

    % 3. Calculate STAs by calling the helper function for each condition
    % --- Stationary STAs ---
    [units(iunit).STA_mean_GREEN_stat, units(iunit).STA_sem_GREEN_stat] = ...
        calculateSTA(stat_SpikeTimes, output.dutyCycle_ch1, output.dutyCycleTime, nSampsBefore);
    
    [units(iunit).STA_mean_UV_stat, units(iunit).STA_sem_UV_stat] = ...
        calculateSTA(stat_SpikeTimes, output.dutyCycle_ch2, output.dutyCycleTime, nSampsBefore);

    % --- Locomotion STAs ---
    [units(iunit).STA_mean_GREEN_run, units(iunit).STA_sem_GREEN_run] = ...
        calculateSTA(loco_SpikeTimes, output.dutyCycle_ch1, output.dutyCycleTime, nSampsBefore);

    [units(iunit).STA_mean_UV_run, units(iunit).STA_sem_UV_run] = ...
        calculateSTA(loco_SpikeTimes, output.dutyCycle_ch2, output.dutyCycleTime, nSampsBefore);
end

toc


%% plot STAs


goodUnits = units([units.rp_contamination]<0.1 & ...
                  [units.amplitude_median]<-50 & ...
                  [units.amplitude_cutoff]<0.1 & ...
                  [units.firing_rate]>=3);

nGood = numel(goodUnits)

for iunit = 1:numel(goodUnits)
    figure, hold on

    sgtitle(['White noise STA for good unit: ', num2str(iunit)])
    xvals = (1:nSampsBefore)./niSampRate;

    subplot(211)
    title('green')
    shadedErrorBar(xvals, goodUnits(iunit).STA_mean_GREEN_stat, goodUnits(iunit).STA_sem_GREEN_stat, 'lineProps', 'k')
    shadedErrorBar(xvals, goodUnits(iunit).STA_mean_GREEN_run, goodUnits(iunit).STA_sem_GREEN_run, 'lineProps', 'r')

    ax = gca;
    ax.XTick = (0:0.02:0.2);
    ax.XTickLabel = -200:20:0;
    ylabel('PWM duty cycle')
    xlabel('Time before spike (ms)')
    xlim([0 0.2])
    defaultAxesProperties(gca, true)

    subplot(212)
    title('UV')
    shadedErrorBar(xvals, goodUnits(iunit).STA_mean_UV_stat, goodUnits(iunit).STA_sem_UV_stat, 'lineProps', 'k')
    shadedErrorBar(xvals, goodUnits(iunit).STA_mean_UV_run, goodUnits(iunit).STA_sem_UV_run, 'lineProps', 'r')

    ax = gca;
    ax.XTick = (0:0.02:0.2);
    ax.XTickLabel = -200:20:0;
    ylabel('PWM duty cycle')
    xlabel('Time before spike (ms)')
    xlim([0 0.2])
    defaultAxesProperties(gca, true)
    pause
    close

end

