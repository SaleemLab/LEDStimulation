function [sta, stc_eigenvectors, stc_eigenvalues] = calculateSTC(spikeTimes, stimulusSignal, stimulusTime, nLags)
% Calculates the Spike-Triggered Average (STA) and performs Spike-Triggered
% Covariance (STC) analysis.
%
% Outputs:
%   sta                - The Spike-Triggered Average vector.
%   stc_eigenvectors   - Matrix whose columns are eigenvectors of the STC matrix.
%   stc_eigenvalues    - Vector of corresponding eigenvalues.

    % --- Handle cases with too few spikes ---
    if numel(spikeTimes) < nLags
        sta = nan(1, nLags);
        stc_eigenvectors = nan(nLags, nLags);
        stc_eigenvalues = nan(nLags, 1);
        return;
    end

    % --- Create Spike-Triggered Ensemble ---
    [~, spike_idx] = histc(spikeTimes, stimulusTime);
    spike_idx(spike_idx <= nLags) = []; % Ensure spikes are not too close to the start
    
    num_spikes = numel(spike_idx);
    ensemble_indices = spike_idx(:) - (nLags - 1) + (0:(nLags - 1));
    spike_triggered_ensemble = stimulusSignal(ensemble_indices);
    
    % --- Calculate STA ---
    sta = mean(spike_triggered_ensemble, 1);
    
    % --- Perform STC Analysis ---
    % Calculate the covariance matrix of the spike-triggered stimulus ensemble.
    stc_matrix = cov(spike_triggered_ensemble);
    
    % Perform eigenvalue decomposition to find principal components.
    [eig_vecs, eig_vals_matrix] = eig(stc_matrix);
    
    % Sort eigenvectors by their corresponding eigenvalues in descending order
    [stc_eigenvalues, sort_idx] = sort(diag(eig_vals_matrix), 'descend');
    stc_eigenvectors = eig_vecs(:, sort_idx);
end