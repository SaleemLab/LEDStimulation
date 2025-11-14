function plot_unit_phase_summary(unit_id, iColor, iLum, unit_phase_analysis, uniqueFreqs)
% PLOT_UNIT_PHASE_SUMMARY Generates a phase-based summary plot for a single neuron.
%
%   Usage:
%   plot_unit_phase_summary(unit_id, iColor, iLum, ...
%                             unit_phase_analysis, uniqueFreqs)
%
%   Inputs:
%   - unit_id:    Index of the unit to plot (e.g., 115)
%   - iColor:     Index for Color (e.g., 1 for UV, 2 for Green)
%   - iLum:       Index for Luminance (e.g., 1 for Dim, 2 for Bright)
%   - unit_phase_analysis: The struct array with phase data
%   - uniqueFreqs: The vector of frequencies for the x-axis (e.g., [2, 5, 10...])
%

% --- 1. Setup Labels and Titles (Customize these to match your indices) ---
try
    colorStrs = {'UV (Contrast -1)', 'Green (Contrast 1)'};
    lumStrs = {'Dim', 'Bright'}; % 1=Dim, 2=Bright
    
    figTitle = sprintf('Unit %d: Phase Analysis | %s | %s (Stat vs Run)', ...
        unit_id, colorStrs{iColor}, lumStrs{iLum});
catch
    % Fallback if indices are out of range
    figTitle = sprintf('Unit %d: Phase Conds: %d, %d (Stat vs Run)', unit_id, iColor, iLum);
end

figure('Name', figTitle, 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.1 0.1 0.7 0.7]);
sgtitle(figTitle, 'FontSize', 14, 'FontWeight', 'bold'); % Super title

% --- 2. Plot 1: Vector Strength (Temporal Fidelity) ---
subplot(2, 2, 1);
try
    % Get data from the new, clean phase struct
    unit_data = unit_phase_analysis(unit_id).conditions;
    
    % Squeeze to get the tuning curve for this condition
    % Dims are [nColors, nLums, nFreqs, nStates]
    
    % --- Stationary Data (iState = 0) ---
    iState_stat = 0;
    % Concatenate struct field with [] before squeezing
    vs_curve_stat = squeeze([unit_data(iColor, iLum, :, iState_stat + 1).vector_strength]);
    
    % --- Running Data (iState = 1) ---
    iState_run = 1;
    % Concatenate struct field with [] before squeezing
    vs_curve_run = squeeze([unit_data(iColor, iLum, :, iState_run + 1).vector_strength]);
    
    % --- Plot ---
    hold on;
    plot(uniqueFreqs, vs_curve_stat, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Stationary');
    plot(uniqueFreqs, vs_curve_run, 'o--', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Running');
    hold off;
    
    title('Vector Strength (Fidelity)');
    xlabel('Frequency (Hz)');
    ylabel('Vector Strength (0-1)');
    grid on;
    legend('show');
    ylim([0 1]); % Vector strength is always between 0 and 1
    xlim([min(uniqueFreqs) max(uniqueFreqs)]);
catch e
    title('Vector Strength (Fidelity)');
    text(0.5, 0.5, 'Error plotting phase data', 'HorizontalAlignment', 'center');
    fprintf('Error plotting Vector Strength: %s\n', e.message);
end

% --- 3. Plot 2: Mean Phase (Latency) ---
subplot(2, 2, 2);
try
    % Get data from the new, clean phase struct
    unit_data = unit_phase_analysis(unit_id).conditions;
    
    % --- Stationary Data (iState = 0) ---
    iState_stat = 0;
    % Concatenate struct field with [] before squeezing
    phase_curve_rad_stat = squeeze([unit_data(iColor, iLum, :, iState_stat + 1).mean_phase_rad]);
    % Use unwrap to prevent 360-degree jumps
    phase_curve_deg_stat = rad2deg(unwrap(phase_curve_rad_stat)); 
    
    % --- Running Data (iState = 1) ---
    iState_run = 1;
    % Concatenate struct field with [] before squeezing
    phase_curve_rad_run = squeeze([unit_data(iColor, iLum, :, iState_run + 1).mean_phase_rad]);
    phase_curve_deg_run = rad2deg(unwrap(phase_curve_rad_run)); 
    
    % --- Plot ---
    hold on;
    plot(uniqueFreqs, phase_curve_deg_stat, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Stationary');
    plot(uniqueFreqs, phase_curve_deg_run, 'o--', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Running');
    hold off;
    
    title('Mean Phase (Latency)');
    xlabel('Frequency (Hz)');
    ylabel('Mean Phase (Degrees)');
    grid on;
    legend('show');
    xlim([min(uniqueFreqs) max(uniqueFreqs)]);
catch e
    title('Mean Phase (Latency)');
    text(0.5, 0.5, 'Error plotting phase data', 'HorizontalAlignment', 'center');
    fprintf('Error plotting Mean Phase: %s\n', e.message);
end

% --- 4. Plot 3: Polar Phase Distribution (Stationary) ---
% *** FIX: Use polaraxes instead of subplot ***
% This creates a polar plot in the bottom-left quadrant
try
    ax_stat = polaraxes('Position', [0.13 0.1 0.33 0.33]);
    iState_stat = 0; % 0 = Stationary
    plot_polar_histograms(ax_stat, unit_phase_analysis(unit_id).conditions, iColor, iLum, iState_stat + 1, uniqueFreqs);
    title('Spike Phase Distributions (Stationary)');
catch e
    title('Spike Phase Distributions (Stationary)');
    text(0.5, 0.5, 'Error plotting polar data', 'HorizontalAlignment', 'center');
    fprintf('Error plotting Polar (Stat): %s\n', e.message);
end


% --- 5. Plot 4: Polar Phase Distribution (Running) ---
% *** FIX: Use polaraxes instead of subplot ***
% This creates a polar plot in the bottom-right quadrant
try
    ax_run = polaraxes('Position', [0.57 0.1 0.33 0.33]);
    iState_run = 1; % 1 = Running
    plot_polar_histograms(ax_run, unit_phase_analysis(unit_id).conditions, iColor, iLum, iState_run + 1, uniqueFreqs);
    title('Spike Phase Distributions (Running)');
catch e
    title('Spike Phase Distributions (Running)');
    text(0.5, 0.5, 'Error plotting polar data', 'HorizontalAlignment', 'center');
    fprintf('Error plotting Polar (Run): %s\n', e.message);
end

end % End of main function


% --- Helper function to plot overlaid polar histograms ---
function plot_polar_histograms(ax, unit_conditions, iColor, iLum, iState_idx, uniqueFreqs)
% This nested function plots all frequency distributions on a single polar axis

% *** FIX: Removed 'axes(ax)' as 'ax' is passed directly to polarhistogram ***
nFreqs = length(uniqueFreqs);
% Use a good colormap
colors = parula(nFreqs); 
legend_handles = gobjects(nFreqs, 1); % Handles for legend
max_r = 0; % To scale all histograms

hold(ax, 'on');

for iFreq = 1:nFreqs
    all_phases = unit_conditions(iColor, iLum, iFreq, iState_idx).all_spike_phases;
    
    % Create legend string
    legend_str = sprintf('%.0f Hz', uniqueFreqs(iFreq));
    
    if ~isempty(all_phases)
        % Create the histogram
        h = polarhistogram(ax, all_phases, 24, 'Normalization', 'probability', ...
                           'FaceColor', colors(iFreq, :), 'FaceAlpha', 0.4, ...
                           'EdgeColor', 'none', 'DisplayName', legend_str);
        
        % Track the max bar height to set axis limits
        max_r = max(max_r, max(h.Values));
        legend_handles(iFreq) = h;
    else
        % Plot a dummy entry just for the legend if no spikes
        legend_handles(iFreq) = polarplot(ax, [0 0], [0 0], 'Color', colors(iFreq,:), ...
                                    'LineWidth', 2, 'DisplayName', [legend_str ' (no spikes)']);
    end
end

% --- Format the polar plot ---
ax.ThetaZeroLocation = 'right'; % Set 0 degrees (phase=0) to the right
ax.ThetaDir = 'counterclockwise';
if max_r == 0, max_r = 1; end % Handle case with no spikes at all
ax.RLim = [0, max_r * 1.15]; % Set radial limit with 15% buffer
legend(ax, legend_handles(isgraphics(legend_handles)), 'Location', 'northeastoutside'); % Only show valid legend entries
hold(ax, 'off');

end % End of helper function

