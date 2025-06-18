%% get session info

Subject = 'M25065';
Session = '20250605';
AcquisitionsID = '0';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);


%% load units
SorterName = 'kilosort4';
units = createClusterTable(ephysDir, SorterName);
units = table2struct(units);


% nidqDir = 'Z:\ibn-vision\DATA\SUBJECTS\M24019\ephys\20240718\nidq_processed';

%% get wheel data from nidq

% nidqDir = 'Z:\ibn-vision\DATA\SUBJECTS\M24077\ephys\20241219\nidq_processed';
% samplingFreq=100; % (resample freq in Hz)
% wheelChans = [3,4]; % NIDQ channels for wheel signal
% wheel = getWheelPos(nidqDir,wheelChans);
% wheel = processNidqWheel(wheel,samplingFreq,'gaussian',0.175); % window size in s

%% get stim timing and async pulse info from nidq
stimONChan = [5];
stimCHANGEChan = [6];
AsyncPulseChan = [1];

nidqMetaFile =  dir(fullfile(nidqDir,'*.meta'));
nidqBinFile =  dir(fullfile(nidqDir,'*.bin'));

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name,nidqBinFile.folder);

stimON = NidqBin(stimONChan,:);
stimCHANGE = NidqBin(stimCHANGEChan,:);
NidqSyncChan = NidqBin(AsyncPulseChan,:);


niSampRate = str2double(NidqMeta.niSampRate);

NidqTime = (1:size(NidqBin,2))./niSampRate;

clear NidqBin

%% get stimulus change times
stimIdx = 1;

if stimIdx<height(ExpInfo)
    nidqidx = [ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx+1)];
else
    nidqidx = [ExpInfo.smp_nidq(stimIdx):numel(NidqTime)];
end


% process nidq data
thisStimON = stimON(nidqidx);
thisStimCHANGE = stimCHANGE(nidqidx);
thisTime = NidqTime(nidqidx);

% get sync pulse times
thisSync = NidqSyncChan(nidqidx);
threshold = mean(thisSync);
thisSync(thisSync<threshold)=0;
thisSync(thisSync>=threshold)=1;
nidq_sync_idx = find(diff(thisSync)==1);
syncTimes_nidq = thisTime(nidq_sync_idx+1);

stimTimeIdx = nidqidx;

% convert to a digital signals
Stim_ON_binary = thisStimON;
threshold = mean(thisStimON);
Stim_ON_binary(Stim_ON_binary<threshold)=0;
Stim_ON_binary(Stim_ON_binary>=threshold)=1;

Stim_CHANGE_binary = thisStimCHANGE;
threshold = mean(thisStimCHANGE);
Stim_CHANGE_binary(Stim_CHANGE_binary<threshold)=0;
Stim_CHANGE_binary(Stim_CHANGE_binary>=threshold)=1;

% find stim on/off (goes from 0 to 1, 1 to 0)
idx = find(diff(Stim_ON_binary)==1); % indexes of where PD goes up
idx = idx+1;
stimONidx = idx;
stimONTime = thisTime(idx); % time stamps of where PD goes up
clear idx

idx = find(diff(Stim_ON_binary)==-1); % indexes of where PD goes up
idx = idx+1;
stimOFFTime = thisTime(idx); % time stamps of where PD goes up
stimOFFidx = idx;
clear idx

% find stim change times (diff == -1 or 1)
idx_up = find(diff(Stim_CHANGE_binary)==1)+1;
idx_down = find(diff(Stim_CHANGE_binary)==-1)+1;
all_idx = sort([idx_up, idx_down, stimOFFidx]);
stimCHANGETimes = thisTime(all_idx);

%% load and process arduino stim file
stimFile = ExpInfo.ArduinoStimFiles{stimIdx};
rawStimTbl = readtable(stimFile);
TOP_idx = find(contains(rawStimTbl.ArduinoString, "TOP"));
TOP_value = str2num(cell2mat(regexp(rawStimTbl.ArduinoString{TOP_idx}, '\d+', 'match')));

C = rawStimTbl.ArduinoString(TOP_idx+1:end-1); % after TOP, ocr values are printed. last element is a flag (-1) so we skip that

% Step 1 & 2: Convert each string to numeric row vector and combine
numCells = cellfun(@(s) str2double(strsplit(s, ',')), C, 'UniformOutput', false);

% Check if all rows have length 3 (if not, error)
rowLengths = cellfun(@numel, numCells);
if any(rowLengths ~= 3)
    error('Not all rows have exactly 3 numbers.');
end

% Convert to numeric matrix (n-by-3)
numMatrix = cell2mat(numCells);

% Step 3: Create table with 3 variables (columns)
StimTable = array2table(numMatrix, 'VariableNames', {'StimIndex', 'Ch1', 'Ch2'});
StimTable.Ch1 = StimTable.Ch1./TOP_value;
StimTable.Ch2 = StimTable.Ch2./TOP_value;

dutyCycle_ch1 = repelem(StimTable.Ch1,diff(all_idx));
dutyCycle_ch2 = repelem(StimTable.Ch2,diff(all_idx));
dutyCycleTime = thisTime(all_idx(1):all_idx(end)-1);

clear stimCHANGETimes all_idx idx_up idx_down Stim_CHANGE_binary Stim_ON_binary

%% load and process wheel data

wheel_tbl = readtable(ExpInfo.WheelFiles{stimIdx});
wheel_tbl.dist = wheel2unit(wheel_tbl.Wheel,4096, 19.5);
wheel_tbl.ddist_dt = [nan; diff(wheel_tbl.dist)./diff(wheel_tbl.ArduinoTime/1000)];
wheel_tbl.ArduinoTime_mid = movmean(wheel_tbl.ArduinoTime/1000,2);
wheel_tbl(1,:)=[];

syncTimes_bonsai = unique(wheel_tbl.LastSyncPulseTime)/1000;

wheel_tbl.nidqTime = mapTimestampsUsingAsyncPulse(syncTimes_bonsai,syncTimes_nidq,wheel_tbl.ArduinoTime_mid);

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

clear numCells numMatrix nidqidx nidqBinFile NidqSyncChan NidqTime rawStimTbl rowLengths stimCHANGE
clear stimON StimTable stimTimeIdx thisStimCHANGE thisStimON thisSync thisTime
clear wheel_zscoreSpeed wheelSpeed wheelTime 


%% generate STAs
tic
timeWindowBefore = 0.2;
bufferVal = 0.5; % extra time buffer
nSampsBefore = ceil(timeWindowBefore * niSampRate);

for iunit = 1:numel(units)
    disp(['Processing Unit ', num2str(iunit)]);

    % Extract spike times and apply filtering criteria
    spikeTimes = units(iunit).spike_times;
    spikeTimes = spikeTimes(spikeTimes > dutyCycleTime(1) + timeWindowBefore + bufferVal);
    spikeTimes = spikeTimes(spikeTimes < dutyCycleTime(end) - timeWindowBefore - bufferVal);

    % sort spikes into stationary and locomotion epochs
    times = spikeTimes(:);
    is_within = (times >= stationary_intervals(:,1)') & (times <= stationary_intervals(:,2)');
    keep_mask = any(is_within, 2);
    kept_times = times(keep_mask);
    stat_SpikeTimes = kept_times;

    times = spikeTimes(:);
    is_within = (times >= locomotion_intervals(:,1)') & (times <= locomotion_intervals(:,2)');
    keep_mask = any(is_within, 2);
    kept_times = times(keep_mask);
    loco_SpikeTimes = kept_times;


    if ~isempty(stat_SpikeTimes)

        numSpikes = numel(stat_SpikeTimes);
        % Find indices of spike times in dutyCycleTime (vectorized search)
        [~, spike_idx] = histc(stat_SpikeTimes, dutyCycleTime);

        % Generate index ranges for all spikes at once (vectorized)
        sample_indices = spike_idx - (nSampsBefore - 1) + (0:(nSampsBefore - 1));

        % Ensure indices are within bounds
        valid_mask = all(sample_indices > 0 & sample_indices <= numel(dutyCycle_ch1), 2);

        STA = nan(numSpikes, nSampsBefore); % Preallocate for speed
        % Extract values in one step (vectorized)
        STA(valid_mask, :) = dutyCycle_ch1(sample_indices(valid_mask, :));

        units(iunit).STA_mean_GREEN_stat = mean(STA,1);
        units(iunit).STA_sem_GREEN_stat = sem(STA,1);

        % Ensure indices are within bounds
        valid_mask = all(sample_indices > 0 & sample_indices <= numel(dutyCycle_ch2), 2);

        STA = nan(numSpikes, nSampsBefore); % Preallocate for speed
        % Extract values in one step (vectorized)
        STA(valid_mask, :) = dutyCycle_ch2(sample_indices(valid_mask, :));

        units(iunit).STA_mean_UV_stat = mean(STA,1);
        units(iunit).STA_sem_UV_stat = sem(STA,1);

    else
        units(iunit).STA_mean_GREEN_stat = nan(1,nSampsBefore);
        units(iunit).STA_sem_GREEN_stat = nan(1,nSampsBefore);
        units(iunit).STA_mean_UV_stat = nan(1,nSampsBefore);
        units(iunit).STA_sem_UV_stat = nan(1,nSampsBefore);
    end

    if ~isempty(loco_SpikeTimes)

        numSpikes = numel(loco_SpikeTimes);
        % Find indices of spike times in dutyCycleTime (vectorized search)
        [~, spike_idx] = histc(loco_SpikeTimes, dutyCycleTime);

        % Generate index ranges for all spikes at once (vectorized)
        sample_indices = spike_idx - (nSampsBefore - 1) + (0:(nSampsBefore - 1));

        % Ensure indices are within bounds
        valid_mask = all(sample_indices > 0 & sample_indices <= numel(dutyCycle_ch1), 2);

        STA = nan(numSpikes, nSampsBefore); % Preallocate for speed
        % Extract values in one step (vectorized)
        STA(valid_mask, :) = dutyCycle_ch1(sample_indices(valid_mask, :));

        units(iunit).STA_mean_GREEN_run = mean(STA,1);
        units(iunit).STA_sem_GREEN_run = sem(STA,1);

        % Ensure indices are within bounds
        valid_mask = all(sample_indices > 0 & sample_indices <= numel(dutyCycle_ch2), 2);

        STA = nan(numSpikes, nSampsBefore); % Preallocate for speed
        % Extract values in one step (vectorized)
        STA(valid_mask, :) = dutyCycle_ch2(sample_indices(valid_mask, :));

        units(iunit).STA_mean_UV_run = mean(STA,1);
        units(iunit).STA_sem_UV_run = sem(STA,1);

    else
        units(iunit).STA_mean_GREEN_run = nan(1,nSampsBefore);
        units(iunit).STA_sem_GREEN_run = nan(1,nSampsBefore);
        units(iunit).STA_mean_UV_run = nan(1,nSampsBefore);
        units(iunit).STA_sem_UV_run = nan(1,nSampsBefore);
    end

end

toc



%% plot STAs

for iunit = 1:numel(units)
    figure, hold on

    title(['White noise STA for unit: ', num2str(iunit)])
    xvals = (1:nSampsBefore)./niSampRate;

    subplot(211)
    title('green')
    shadedErrorBar(xvals, units(iunit).STA_mean_GREEN_stat, units(iunit).STA_sem_GREEN_stat, 'lineProps', 'k')
    shadedErrorBar(xvals, units(iunit).STA_mean_GREEN_run, units(iunit).STA_sem_GREEN_run, 'lineProps', 'r')
    
    ax = gca;
    ax.XTick = (0:0.02:0.2);
    ax.XTickLabel = -200:20:0;
    ylabel('PWM duty cycle')
    xlabel('Time before spike (ms)')
    xlim([0 0.2])
    defaultAxesProperties(gca, true)

    subplot(212)
    title('UV')
    shadedErrorBar(xvals, units(iunit).STA_mean_UV_stat, units(iunit).STA_sem_UV_stat, 'lineProps', 'k')
    shadedErrorBar(xvals, units(iunit).STA_mean_UV_run, units(iunit).STA_sem_UV_run, 'lineProps', 'r')

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

