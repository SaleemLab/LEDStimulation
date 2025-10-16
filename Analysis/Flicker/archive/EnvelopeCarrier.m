%% get session info

Subject = 'M24077';
Session = '20241219';
AcquisitionsID = '2';
BaseDir = 'Z:\ibn-vision\DATA\SUBJECTS';
% kilosortDir = fullfile(BaseDir,Subject,'ephys',Session,'spike_sorting','probe0','sorters','kilosort3_merged');

[ExpInfo, ephysDir, nidqDir, bonsaiDir] = getExpInfo(BaseDir, Subject, Session,AcquisitionsID);


%% get wheel data from nidq
samplingFreq=100; % (resample freq in Hz)
wheelChans = [3,4]; % NIDQ channels for wheel signal
wheel = getWheelPos(nidqDir,wheelChans);
wheel = processNidqWheel(wheel,samplingFreq,'gaussian',0.175); % window size in s

%% get stim on times
stimONChan = [5];

nidqMetaFile =  dir(fullfile(nidqDir,'*.meta'));
nidqBinFile =  dir(fullfile(nidqDir,'*.bin'));

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(1, inf, NidqMeta, nidqBinFile.name,nidqBinFile.folder);

stimON = NidqBin(stimONChan,:);

niSampRate = str2double(NidqMeta.niSampRate);

NidqTime = (1:size(NidqBin,2))./niSampRate;

clear NidqBin

%% load units
SorterName = 'kilosort3_merged';
units = createClusterTable(ephysDir, SorterName);
units = table2struct(units);


%% get stimulus start times for flicker

stimsToCompare = {10} %{3, 4, 11, 12}; % each cell is vector of stim idx from ExpInfo

for iStimToCompare = 1:numel(stimsToCompare)

stim_idx_to_use = stimsToCompare{iStimToCompare};

ExpInfo.Stimulus(stim_idx_to_use(1))
allStimTimes = [];
allStimTimesImec0=[];

for stim_index = stim_idx_to_use

    stimIdx = stim_index;

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
    idx2 = find(diffDiffTimes>5); % For new trials its ~2, i.e. >1.5s
    idx2=idx2+1;
    trialStartIdx = idx(idx2); % convert back to original PD index
    trialStartIdx = [idx(1), trialStartIdx]; % add the first trial
    stimTimes = thisTime(trialStartIdx); % get the PD timestamp for trial start idx

    % timeFrame = [NidqTime(stimTimeIdx(1)), NidqTime(stimTimeIdx(end))-5];
    % 
    % % convert stim times to imec0 timestamp using interp1
    % % (maybe this should be converted to tprime)
    % this_nidqsync_nidq = nidq_syncTimes_nidq(nidq_syncTimes_nidq >= timeFrame(1) & nidq_syncTimes_nidq <= timeFrame(2));
    % this_nidqsync_imec0 = nidq_syncTimes_imec0(nidq_syncTimes_imec0 >= timeFrame(1) & nidq_syncTimes_imec0 <= timeFrame(2));
    % this_nidqsync_nidq=this_nidqsync_nidq(:); this_nidqsync_imec0=this_nidqsync_imec0(:);
    % 
    % 
    % while numel(this_nidqsync_imec0)>numel(this_nidqsync_nidq)
    %     dv = this_nidqsync_imec0(1:numel(this_nidqsync_nidq))-this_nidqsync_nidq;
    %     idx = find(abs(dv)>0.5,1,'first');
    %     this_nidqsync_imec0(idx)=[];
    % end
    % 
    % while numel(this_nidqsync_nidq)>numel(this_nidqsync_imec0)
    %     dv = this_nidqsync_nidq(1:numel(this_nidqsync_imec0))-this_nidqsync_imec0;
    %     idx = find(abs(dv)>0.5,1,'first');
    %     this_nidqsync_nidq(idx)=[];
    % end
    % 
    % 
    % stimTimes_imec0 = interp1(this_nidqsync_nidq, this_nidqsync_imec0, stimTimes);

    allStimTimes = cat(2,allStimTimes, stimTimes);
    % allStimTimesImec0 = cat(2,allStimTimesImec0, stimTimes_imec0);


end


% get trial params
% quick hack for photodidoe
photodiode.up_times = allStimTimes;
photodiode.down_times = photodiode.up_times+6;

allFlickerTrials = createTrialsStruct(ExpInfo,stim_idx_to_use,photodiode,wheel);

stat_idx = [allFlickerTrials.meanWheelSpeed]<=1;
run_idx = [allFlickerTrials.meanWheelSpeed]>3;
[allFlickerTrials.runFlag] = deal([nan]);
[allFlickerTrials(stat_idx).runFlag] = deal([0]);
[allFlickerTrials(run_idx).runFlag] = deal([1]);


% get set of intervals for trials
statTrials = allFlickerTrials([allFlickerTrials.runFlag]==0);
runTrials = allFlickerTrials([allFlickerTrials.runFlag]==1);

flickFreqs = unique([allFlickerTrials.CarrierFrequency]);

for ifreq = 1:numel(flickFreqs)
    % stationary/locomotion trials
    stimTrials = statTrials([statTrials.CarrierFrequency]==flickFreqs(ifreq));
    stat_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;

    % stimTrials = runTrials([runTrials.CarrierFrequency]==flickFreqs(ifreq));
    % run_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;
end


end
% end

%%

for iunit = 1:numel(units)
    iunit
options=struct;
options.smoothWidth = 250;
options.preTime = 1000;
options.postTime=1000;
figure
[psth, blank, binMidPoints] = makePSTH(units(iunit).spike_times*1000,stat_intervalsCell,[],options);
legend({' ','10', '40', '200'})
title(iunit)
pause
close

end




