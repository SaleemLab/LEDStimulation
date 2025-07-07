function [sta_mean, sta_sem] = calculateSTA(spikeTimes, stimulusSignal, stimulusTime, nSampsBefore)
% Calculates the Spike-Triggered Average (STA) for a given set of spikes and a stimulus signal.
%
% Inputs:
%   spikeTimes      - Vector of spike times.
%   stimulusSignal  - The continuous signal to average (e.g., dutyCycle_ch1).
%   stimulusTime    - The time vector corresponding to the stimulusSignal.
%   nSampsBefore    - The number of samples to include in the window before each spike.
%
% Outputs:
%   sta_mean        - The mean spike-triggered average.
%   sta_sem         - The standard error of the mean of the STA.

    if isempty(spikeTimes)
        sta_mean = nan(1, nSampsBefore);
        sta_sem = nan(1, nSampsBefore);
        return;
    end

    numSpikes = numel(spikeTimes);

    % Find indices of spike times in the stimulus time vector.
    % Note: 'discretize' is the modern replacement for 'histc' but 'histc' is used here for legacy compatibility.
    [~, spike_idx] = histc(spikeTimes, stimulusTime);

    % Generate index ranges for all spikes at once (vectorized approach)
    % This creates a [numSpikes x nSampsBefore] matrix of indices.
    sample_indices = spike_idx(:) - (nSampsBefore - 1) + (0:(nSampsBefore - 1));

    % Create a mask to exclude any windows that go out of bounds (e.g., for spikes too close to the start)
    valid_mask = all(sample_indices > 0 & sample_indices <= numel(stimulusSignal), 2);

    % Pre-allocate a matrix to hold all the stimulus snippets
    STA_matrix = nan(numSpikes, nSampsBefore);

    % Extract stimulus snippets for all valid spikes in a single, vectorized step
    STA_matrix(valid_mask, :) = stimulusSignal(sample_indices(valid_mask, :));

    % Calculate the mean and SEM, ignoring any rows that remained NaN
    sta_mean = mean(STA_matrix, 1, 'omitnan');
    sta_sem = std(STA_matrix, 0, 1, 'omitnan') / sqrt(sum(valid_mask));
end