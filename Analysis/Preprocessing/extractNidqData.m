function nidq = extractNidqData(ExpInfo, igkey, nidqDir, nidq_channels, nidq_1hz_file, nidq_tprime_file, doPlot)
%extractNidqData Extracts, thresholds, and analyzes NIDQ data for a specific stimulus.
%
%   This function reads a segment of a NIDQ binary file corresponding to a
%   given experimental key. It dynamically extracts specified channels,
%   calculates a binary threshold, maps the time base, and finds the times
%   of rising edges on the asynchronous pulse channel.
%
%   SYNTAX:
%   nidq = extractNidqData(ExpInfo, igkey, nidqDir, nidq_channels)
%   nidq = extractNidqData(..., doPlot)
%
%   INPUTS:
%   ExpInfo         - Structure containing experiment information.
%   igkey           - Scalar integer specifying which stimulus to extract.
%   nidqDir         - String path to the directory with NIDQ .meta and .bin files.
%   nidq_channels   - Structure where each field name is a desired output name
%                     (e.g., 'asyncPulse') and its value is the channel number.
%   nidq_1hz_file   - (Optional) Path to NIDQ pulse times file for time mapping.
%   nidq_tprime_file - (Optional) Path to imec0 reference times file for time mapping.
%   doPlot          - (Optional) Logical flag (true/false) to plot channels. Defaults to false.
%
%   OUTPUT:
%   nidq            - Structure containing the extracted data, time vectors,
%                     and the calculated asynchronous pulse times.

%% 1. Argument Handling
if nargin < 7
    doPlot = false;
end

%% 2. Identify Data Segment
% Find the index corresponding to the chosen stimulus key
g_idx = find(ExpInfo.gKey == igkey);
if isempty(g_idx)
    error('Specified igkey (%d) not found in ExpInfo.gKey.', igkey);
end

% Determine the start sample and number of samples to read
samp0 = ExpInfo.smp_nidq(g_idx);
if g_idx ~= height(ExpInfo)
    nSamp = ExpInfo.smp_nidq(g_idx + 1) - samp0;
else
    % If it's the last stimulus, read to the end of the file
    nSamp = Inf;
end

%% 3. Read Raw Data from Binary File
nidqMetaFile = dir(fullfile(nidqDir, '*.meta'));
nidqBinFile = dir(fullfile(nidqDir, '*.bin'));

if isempty(nidqMetaFile) || isempty(nidqBinFile)
    error('Could not find .meta or .bin file in the specified directory: %s', nidqDir);
end

NidqMeta = ReadMeta(fullfile(nidqMetaFile.folder, nidqMetaFile.name));
NidqBin = ReadBin(samp0, nSamp, NidqMeta, nidqBinFile.name, nidqBinFile.folder);

%% 4. Extract and Threshold Channels
channel_fields = fieldnames(nidq_channels);
if isempty(channel_fields)
    error('The nidq_channels struct cannot be empty.');
end

for i = 1:length(channel_fields)
    fieldName = channel_fields{i};
    channel_idx = nidq_channels.(fieldName);

    rawData = NidqBin(channel_idx, :);

    % Automatically determine and apply a binary threshold
    if range(rawData) < eps % Handle flat signals
        nidq.(fieldName) = zeros(size(rawData), 'double');
    else
        % Calculate threshold as the midpoint between the 5th and 95th percentiles.
        quantiles = quantile(double(rawData), [0.05, 0.95]);
        threshold = mean(quantiles);
        nidq.(fieldName) = double(rawData > threshold);
    end
end

clear NidqBin;

%% 5. Generate NIDQ Time Vector
first_channel_name = channel_fields{1};
nSampActual = size(nidq.(first_channel_name), 2);

t0 = ExpInfo.sec_nidq(g_idx);
niSampRate = str2double(NidqMeta.niSampRate);

nidq.time = t0 + (0:nSampActual-1) / niSampRate;

%% 6. Map NIDQ Time to imec0 Time Base (Optional)
if nargin >= 6 && ~isempty(nidq_1hz_file) && ~isempty(nidq_tprime_file)
    t_pulse_sig = readmatrix(nidq_1hz_file);
    t_pulse_ref = readmatrix(nidq_tprime_file);
    nidq.imecTime = interp1(t_pulse_sig, t_pulse_ref, nidq.time, 'linear', 'extrap');
else
    nidq.imecTime = []; % Explicitly set to empty if not mapped
end

%% 7. Plot Channels (Optional)
if doPlot
    % Determine which time vector and label to use for the x-axis
    if ~isempty(nidq.imecTime)
        timeVector = nidq.imecTime;
        xLabelText = 'Time (imec0 base, s)';
    else
        timeVector = nidq.time;
        xLabelText = 'Time (NIDQ base, s)';
    end

    figure;
    num_channels = length(channel_fields);
    ax = gobjects(num_channels, 1); % Pre-allocate axes handles

    for i = 1:num_channels
        fieldName = channel_fields{i};
        ax(i) = subplot(num_channels, 1, i);
        plot(timeVector, nidq.(fieldName));
        
        % Set title, preventing underscores from being interpreted as subscripts
        title(fieldName, 'Interpreter', 'none');
        ylabel('State');
        ylim([-0.1, 1.1]); % Set y-axis limits for binary data
    end
    
    % Label the x-axis on the last subplot and link all axes
    xlabel(xLabelText);
    linkaxes(ax, 'x');
end

%% 8. Find Asynchronous Pulse Times
if isfield(nidq, 'asyncPulse')
    % Find indices where the signal goes from 0 to 1
    rising_edge_indices = find(diff(nidq.asyncPulse) == 1) + 1;

    % Get the timestamps for these events using the NIDQ time base
    nidq.asyncPulseTimes_nidq = nidq.time(rising_edge_indices);

    % Get timestamps using the imec0 time base, if it exists
    if ~isempty(nidq.imecTime)
        nidq.asyncPulseTimes_imec = nidq.imecTime(rising_edge_indices);
    else
        nidq.asyncPulseTimes_imec = [];
    end
else
    % If asyncPulse wasn't requested, create empty fields for consistency
    nidq.asyncPulseTimes_nidq = [];
    nidq.asyncPulseTimes_imec = [];
end

end