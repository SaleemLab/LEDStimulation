%% Spike-Triggered Analysis Suite (STA, STC, ICA)
% This script demonstrates how to perform spike-triggered average (STA),
% spike-triggered covariance (STC), and spike-triggered independent 
% component analysis (ICA) on a simulated neuron responding to a 1D 
% Gaussian white noise stimulus.

% Clear workspace and close all figures
clear;
clc;
close all;

%% 1. SETUP AND DATA SIMULATION
% In a real scenario, you would load your stimulus and spike times here.
% For this demonstration, we simulate data from a Linear-Nonlinear-Poisson
% (LNP) model neuron.

% --- Simulation Parameters ---
fs = 1000;                  % Sampling rate (Hz)
duration_s = 200;           % Duration of the experiment (seconds)
n_samples = duration_s * fs;% Total number of samples
t = (0:n_samples-1) / fs;   % Time vector

% --- Stimulus: Gaussian White Noise ---
stimulus = randn(1, n_samples);
stimulus = (stimulus - mean(stimulus)) / std(stimulus); % Normalize to mean 0, std 1

% --- LNP Model Neuron Simulation ---
% We create a "ground truth" filter that our analysis should ideally recover.
% A Gabor-like filter is a common choice.
filter_time_ms = 150; % How far back in time the filter extends
filter_n_samples = round(filter_time_ms / 1000 * fs);
filter_t = (-filter_n_samples+1:0) / fs;

% The "true" filter of our model neuron
true_filter = exp(-((filter_t*fs)/20).^2) .* cos(2*pi*30*filter_t);
true_filter = true_filter / norm(true_filter); % Normalize

% Linear Step: Convolve stimulus with the filter
generator_signal = conv(stimulus, fliplr(true_filter), 'same');

% Nonlinear Step: Apply a nonlinearity (e.g., squaring)
% This makes the neuron respond to both positive and negative filter outputs.
nonlinearity_output = generator_signal.^2;

% Poisson Firing: Convert to an instantaneous firing rate and generate spikes
baseline_rate = 5;      % Spontaneous firing rate (Hz)
max_driven_rate = 100;  % Max rate driven by stimulus (Hz)
firing_rate = baseline_rate + max_driven_rate * nonlinearity_output;
spikes = rand(1, n_samples) < (firing_rate / fs);
spike_indices = find(spikes);
n_spikes = length(spike_indices);

fprintf('Data generated. Number of spikes: %d\n', n_spikes);

% --- Plot Simulated Data ---
figure('Name', 'Simulated Data Overview');
subplot(3, 1, 1);
plot(t, stimulus);
title('Stimulus (Gaussian White Noise)');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 5]);

subplot(3, 1, 2);
plot(t, firing_rate);
title('Instantaneous Firing Rate');
xlabel('Time (s)'); ylabel('Rate (Hz)');
xlim([0 5]);

subplot(3, 1, 3);
plot(t, spikes, '.');
hold on;
plot(spike_indices/fs, ones(1, n_spikes), 'r|');
title('Generated Spike Train');
xlabel('Time (s)');
xlim([0 5]);
ylim([0 1.1]);


%% 2. SPIKE-TRIGGERED AVERAGE (STA)
% The STA is the average stimulus segment preceding a spike. It's an
% estimate of the neuron's linear receptive field.

fprintf('\nCalculating Spike-Triggered Average (STA)...\n');

% --- Define the time window for analysis ---
time_window_ms = 200; % ms
time_window_samples = round(time_window_ms / 1000 * fs);

% --- Collect the spike-triggered ensemble (STE) ---
% The STE is a matrix where each row is the stimulus segment before a spike.
spike_triggered_ensemble = zeros(n_spikes, time_window_samples);
spike_count = 0;

for i = 1:n_spikes
    spike_idx = spike_indices(i);
    % Ensure the window does not go past the beginning of the stimulus
    if spike_idx > time_window_samples
        spike_count = spike_count + 1;
        stim_segment = stimulus(spike_idx - time_window_samples : spike_idx - 1);
        spike_triggered_ensemble(spike_count, :) = stim_segment;
    end
end
% Trim unused rows if any spikes were too early
spike_triggered_ensemble = spike_triggered_ensemble(1:spike_count, :);

% --- Calculate STA ---
sta = mean(spike_triggered_ensemble, 1);

% --- Plot STA ---
figure('Name', 'Spike-Triggered Average (STA)');
t_window = (-time_window_samples+1:0) / fs * 1000; % Time vector in ms
plot(t_window, sta, 'b-', 'LineWidth', 2, 'DisplayName', 'Calculated STA');
hold on;
% Plot the ground truth filter for comparison (scaled to match STA amplitude)
plot(filter_t*1000, true_filter * (max(abs(sta))/max(abs(true_filter))), 'r--', 'LineWidth', 2, 'DisplayName', 'True Filter (scaled)');
title('Spike-Triggered Average (STA)');
xlabel('Time before spike (ms)');
ylabel('Stimulus Amplitude');
grid on;
legend show;
ax = gca;
ax.XAxis.Direction = 'reverse'; % Time runs backwards from spike at t=0
xlim([min(t_window), 0]);


%% 3. SPIKE-TRIGGERED COVARIANCE (STC)
% STC analyzes the variance structure of the spike-triggered ensemble.
% Eigenvectors of the STC matrix reveal additional stimulus dimensions
% that modulate the neuron's firing rate.

fprintf('Calculating Spike-Triggered Covariance (STC)...\n');

% --- Calculate the STC matrix ---
% Note: We use the same spike_triggered_ensemble from the STA calculation.
stc = cov(spike_triggered_ensemble);

% --- Perform Eigenvalue Decomposition ---
[eigenvectors, eigenvalues_matrix] = eig(stc);
eigenvalues = diag(eigenvalues_matrix);

% Sort eigenvalues and corresponding eigenvectors in descending order
[eigenvalues, sort_idx] = sort(eigenvalues, 'descend');
eigenvectors = eigenvectors(:, sort_idx);

% --- Plot STC Results ---
figure('Name', 'Spike-Triggered Covariance (STC) Analysis');

% Plot the covariance matrix
subplot(2, 2, 1);
imagesc(t_window, t_window, stc);
title('Spike-Triggered Covariance Matrix');
xlabel('Time before spike (ms)');
ylabel('Time before spike (ms)');
colorbar;

% Plot the eigenvalue spectrum
subplot(2, 2, 2);
plot(eigenvalues, 'o-', 'LineWidth', 1.5);
title('Eigenvalue Spectrum');
xlabel('Eigenvalue Rank');
ylabel('Eigenvalue');
grid on;
set(gca, 'YScale', 'log'); % Often better to view on a log scale

% Plot the most significant eigenvectors
% The eigenvector with the largest eigenvalue
subplot(2, 2, 3);
plot(t_window, eigenvectors(:, 1), 'm-', 'LineWidth', 2, 'DisplayName', 'EV 1 (Largest \lambda)');
hold on;
% For our LNP model with a squaring nonlinearity, the STA is near zero,
% but the first eigenvector of the STC should recover the true filter.
plot(filter_t*1000, true_filter * (max(abs(eigenvectors(:, 1)))/max(abs(true_filter))), 'r--', 'LineWidth', 2, 'DisplayName', 'True Filter (scaled)');
title('Most Significant Excitatory Filter');
xlabel('Time before spike (ms)');
ylabel('Amplitude');
ax = gca; ax.XAxis.Direction = 'reverse'; xlim([min(t_window), 0]);
grid on;
legend show;

% The eigenvector with the smallest eigenvalue (most suppressive direction)
subplot(2, 2, 4);
plot(t_window, eigenvectors(:, end), 'g-', 'LineWidth', 2, 'DisplayName', ['EV ' num2str(length(eigenvalues)) ' (Smallest \lambda)']);
title('Most Significant Suppressive Filter');
xlabel('Time before spike (ms)');
ylabel('Amplitude');
ax = gca; ax.XAxis.Direction = 'reverse'; xlim([min(t_window), 0]);
grid on;
legend show;


%% 4. SPIKE-TRIGGERED INDEPENDENT COMPONENT ANALYSIS (ICA)
% ICA attempts to find a set of filters such that their outputs, when
% projected onto the spike-triggered ensemble, are as statistically
% independent as possible. This can reveal non-orthogonal features.
%
% NOTE: This requires the Statistics and Machine Learning Toolbox for 'rica'.

fprintf('\nCalculating Spike-Triggered Independent Components (ICA)...\n');

if license('test', 'statistics_toolbox')
    num_components = 4; % How many independent components to find
    
    % rica (Reconstruction ICA) is a good choice. It finds an orthogonal
    % rotation that maximizes the non-Gaussianity of the components.
    mdl = rica(spike_triggered_ensemble, num_components);
    
    % The learned filters are often stored as columns in the transposed 
    % TransformWeights matrix.
    ica_filters = mdl.TransformWeights; % Transpose the weights matrix
    
    % --- Plot ICA Filters ---
    figure('Name', 'Spike-Triggered Independent Components (ICA)');
    sgtitle('ICA Filters'); % Super title for the figure
    for i = 1:num_components
        subplot(num_components, 1, i);
        % Plot the i-th column, which corresponds to one filter
        plot(t_window, ica_filters(:, i), 'LineWidth', 2);
        title(['Independent Component ' num2str(i)]);
        ylabel('Weight');
        ax = gca; ax.XAxis.Direction = 'reverse'; xlim([min(t_window), 0]);
        grid on;
        if i == num_components
            xlabel('Time before spike (ms)');
        end
    end
else
    fprintf('Statistics and Machine Learning Toolbox not found. Skipping ICA.\n');
end

fprintf('\nAnalysis complete.\n');
