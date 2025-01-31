%% get session info

Subject = 'M24077';
Session = '20241219';
AcquisitionsID = '2';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);


%% load units
SorterName = 'kilosort3_merged';
units = createClusterTable(ephysDir, SorterName);
units = table2struct(units);


% nidqDir = 'Z:\ibn-vision\DATA\SUBJECTS\M24019\ephys\20240718\nidq_processed';

%% get wheel data from nidq

% nidqDir = 'Z:\ibn-vision\DATA\SUBJECTS\M24077\ephys\20241219\nidq_processed';
samplingFreq=100; % (resample freq in Hz)
wheelChans = [3,4]; % NIDQ channels for wheel signal
wheel = getWheelPos(nidqDir,wheelChans);
wheel = processNidqWheel(wheel,samplingFreq,'gaussian',0.175); % window size in s

%% get stim timing info from nidq
stimONChan = [5];
stimCHANGEChan = [6];

nidqMetaFile =  dir(fullfile(nidqDir,'*.meta'));
nidqBinFile =  dir(fullfile(nidqDir,'*.bin'));

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name,nidqBinFile.folder);

stimON = NidqBin(stimONChan,:);
stimCHANGE = NidqBin(stimCHANGEChan,:);

niSampRate = str2double(NidqMeta.niSampRate);

NidqTime = (1:size(NidqBin,2))./niSampRate;

clear NidqBin

%% get stimulus change times
stimIdx = 6;

if stimIdx<height(ExpInfo)
    nidqidx = [ExpInfo.smp_nidq(stimIdx):ExpInfo.smp_nidq(stimIdx+1)];
else
    nidqidx = [ExpInfo.smp_nidq(stimIdx):numel(NidqTime)];
end


% get stimON and stimCHANGE times
thisStimON = stimON(nidqidx);
thisStimCHANGE = stimCHANGE(nidqidx);
thisTime = NidqTime(nidqidx);
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

%% load duty cycle values for each stim change
opts = delimitedTextImportOptions("NumVariables", 3);
% Specify range and delimiter
opts.Delimiter = ",";
opts.VariableTypes = ["string", "double", "double"];

% ardStimLog = readtable(fullfile('D:\Code\LEDStimulation\Analysis\WhiteNoise','M24077_WhiteNoise_UV_ArduinoStimLog2024-12-19T17_51_48'),opts);
ardStimLog = readtable(fullfile('D:\Code\LEDStimulation\Analysis\WhiteNoise','M24077_WhiteNoise_GREEN_ArduinoStimLog2024-12-19T18_23_40'),opts);
TOP_idx = find(contains(ardStimLog.Var1, "TOP"));
TOP_value = double(regexp(ardStimLog.Var1(TOP_idx), '\d+', 'match'));

reqLumVals = double(ardStimLog.Var1(TOP_idx+1:end)); % last value is -1 to show stimulus finished...
if reqLumVals(end)==-1, reqLumVals(end)=[]; end

dutyCycleValues = reqLumVals./TOP_value;

dutyCycleRecon = repelem(dutyCycleValues,diff(all_idx));
dutyCycleTime = thisTime(all_idx(1):all_idx(end)-1);

clear stimCHANGETimes all_idx idx_up idx_down Stim_CHANGE_binary Stim_ON_binary

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

    if ~isempty(spikeTimes)
    numSpikes = numel(spikeTimes);
    STA = nan(numSpikes, nSampsBefore); % Preallocate for speed

    % Find indices of spike times in dutyCycleTime (vectorized search)
    [~, spike_idx] = histc(spikeTimes, dutyCycleTime);

    % Generate index ranges for all spikes at once (vectorized)
    sample_indices = spike_idx - (nSampsBefore - 1) + (0:(nSampsBefore - 1));

    % Ensure indices are within bounds
    valid_mask = all(sample_indices > 0 & sample_indices <= numel(dutyCycleRecon), 2);
    
    % Extract values in one step (vectorized)
    STA(valid_mask, :) = dutyCycleRecon(sample_indices(valid_mask, :));

    units(iunit).STA_mean_GREEN = mean(STA,1);
    units(iunit).STA_sem_GREEN = sem(STA,1);

    else


    units(iunit).STA_mean_GREEN = nan(1,nSampsBefore);
    units(iunit).STA_sem_GREEN = nan(1,nSampsBefore);


    end

end

toc



%% plot STAs

for iunit = 1:numel(units)
    figure, hold on

title(['White noise STA for unit: ', num2str(iunit)])
xvals = (1:nSampsBefore)./niSampRate;

shadedErrorBar(xvals, units(iunit).STA_mean_UV, units(iunit).STA_sem_UV, 'lineProps', 'b')
shadedErrorBar(xvals, units(iunit).STA_mean_GREEN, units(iunit).STA_sem_GREEN, 'lineProps', 'g')

ax = gca;
ax.XTick = (0:0.05:0.2);
ax.XTickLabel = -200:50:0;
ylabel('PWM duty cycle')
xlabel('Time before spike (ms)')
xlim([0 0.2])
defaultAxesProperties(gca, true)
pause
close

end

