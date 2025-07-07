function output = processStimulusData(stimIdx, ExpInfo, nidqDir)
% processStimulusData Loads and processes NIDQ/Arduino data for a specific stimulus.
%
% This function is a complete, memory-efficient pipeline. It locates NIDQ
% files, reads ONLY the required data segment from the binary file on disk,
% extracts event times, and correlates them with parameters from an Arduino log file.
%
% Syntax:
%   output = processStimulusData(stimIdx, ExpInfo, nidqDir)
%
% Inputs:
%   stimIdx         - Scalar index of the stimulus trial to process.
%   ExpInfo         - Table or struct with experiment metadata, including:
%                     .smp_nidq: Start sample index for each stimulus.
%                     .ArduinoStimFiles: Cell array of paths to Arduino log files.
%   nidqDir         - String path to the directory containing the NIDQ
%                     .bin and .meta files for the recording session.
%
% Outputs:
%   output - A structure containing the processed data.

%% 1. Define Constants and Locate NIDQ Files
STIM_ON_CHAN = 5;
STIM_CHANGE_CHAN = 6;
SYNC_PULSE_CHAN = 1;

nidqMetaFile = dir(fullfile(nidqDir, '*.meta'));
nidqBinFile = dir(fullfile(nidqDir, '*.bin'));
if isempty(nidqMetaFile) || isempty(nidqBinFile)
    error('Could not find the required .meta and .bin files in: %s', nidqDir);
end
NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
niSampRate = str2double(NidqMeta.niSampRate);

%% 2. Read Required Data Segment from Disk
start_samp = ExpInfo.smp_nidq(stimIdx);
if stimIdx < height(ExpInfo)
    end_samp_inclusive = ExpInfo.smp_nidq(stimIdx + 1);
    n_samp = end_samp_inclusive - start_samp + 1;
else
    n_samp = inf;
end
segmentData = ReadBin(start_samp, n_samp, NidqMeta, nidqBinFile.name, nidqDir);
actual_n_samp = size(segmentData, 2);
end_samp = start_samp + actual_n_samp;
stimON_segment = segmentData(STIM_ON_CHAN, :);
stimCHANGE_segment = segmentData(STIM_CHANGE_CHAN, :);
sync_segment = segmentData(SYNC_PULSE_CHAN, :);
clear segmentData;
time_indices = start_samp:(end_samp-1);
thisTime = time_indices ./ niSampRate;
output.stimTimeIdx = time_indices;

%% 3. Process NIDQ Signals to Find Event Times
% Get sync pulse times
is_high = sync_segment >= mean(sync_segment);
nidq_sync_idx = find(diff(is_high) == 1) + 1;
output.syncTimes_nidq = thisTime(nidq_sync_idx);

% Get stimulus ON/OFF times
is_on = stimON_segment >= mean(stimON_segment);
on_idx_local = find(diff(is_on) == 1) + 1;
output.stimONTime = thisTime(on_idx_local);
off_idx_local = find(diff(is_on) == -1) + 1;
output.stimOFFTime = thisTime(off_idx_local);

% Get stimulus CHANGE times
is_change = stimCHANGE_segment >= mean(stimCHANGE_segment);
idx_up_local = find(diff(is_change) == 1) + 1;
idx_down_local = find(diff(is_change) == -1) + 1;

% The stimCHANGE channel is forced low when the stimON channel goes low,
% marking the end of a stimulus presentation. This particular stimCHANGE
% falling edge is an artifact and should not be counted as a stimulus event.
if ~isempty(off_idx_local)
    % Identify the index of the final stimulus OFF event within this segment.
    final_stim_off_idx = max(off_idx_local);
    
    % Find if a stimCHANGE falling edge occurred at the exact same index.
    is_artifact = (idx_down_local == final_stim_off_idx);
    
    % Remove that specific artifactual event from the list of change events.
    idx_down_local(is_artifact) = [];
end

all_idx_local = sort([idx_up_local, idx_down_local, off_idx_local]);

%% 4. Load and Process Arduino Stimulus File
stimFile = ExpInfo.ArduinoStimFiles{stimIdx};
rawStimTbl = readtable(stimFile);
top_row_idx = find(contains(rawStimTbl.ArduinoString, "TOP"), 1);
if isempty(top_row_idx)
    error('The "TOP" marker was not found in the Arduino file: %s', stimFile);
end
top_value_str = regexp(rawStimTbl.ArduinoString{top_row_idx}, '\d+', 'match');
top_value = str2double(cell2mat(top_value_str));
data_strings = rawStimTbl.ArduinoString(top_row_idx+1 : end-1);
numRows = numel(data_strings);
numMatrix = zeros(numRows, 3);
for i = 1:numRows
    parsed_row = sscanf(data_strings{i}, '%f,%f,%f');
    if numel(parsed_row) == 3
        numMatrix(i, :) = parsed_row';
    end
end
output.StimTable = array2table(numMatrix, 'VariableNames', {'StimIndex', 'Ch1', 'Ch2'});
output.StimTable.Ch1 = output.StimTable.Ch1 / top_value;
output.StimTable.Ch2 = output.StimTable.Ch2 / top_value;

%% 5. Generate Reconstructed Duty Cycle Signals
if numel(all_idx_local) ~= (height(output.StimTable) + 1)
     warning('Mismatch between NIDQ events (%d) and Arduino stimuli (%d). Timing may be inaccurate.', numel(all_idx_local), height(output.StimTable));
end
repeats = diff(all_idx_local);
output.dutyCycle_ch1 = repelem(output.StimTable.Ch1, repeats);
output.dutyCycle_ch2 = repelem(output.StimTable.Ch2, repeats);
start_idx = all_idx_local(1);
end_idx = start_idx + sum(repeats) - 1;
output.dutyCycleTime = thisTime(start_idx:end_idx);

end