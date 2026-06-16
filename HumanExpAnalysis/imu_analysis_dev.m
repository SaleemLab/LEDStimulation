% 1. Create a clean time axis in seconds
t = (imuData.TimeUs - imuData.TimeUs(1)) / 1e6;

% ==========================================
% PART 1: ORIENTATION (YAW & ROLL)
% ==========================================
q_w = imuData.QuatReal;
q_x = imuData.QuatI;
q_y = imuData.QuatJ;
q_z = imuData.QuatK;

% Yaw (Heading)
heading_deg = atan2(2 .* (q_w .* q_z + q_x .* q_y), 1 - 2 .* (q_y.^2 + q_z.^2)) .* (180/pi);
wrap_indices = find(abs(diff(heading_deg)) > 300);
heading_plot = heading_deg;
heading_plot(wrap_indices) = NaN;

% Roll (Side-to-Side Sway) centered at 0
roll_deg = atan2(2 .* (q_w .* q_x + q_y .* q_z), 1 - 2 .* (q_x.^2 + q_y.^2)) .* (180/pi);
roll_sway = roll_deg - mean(roll_deg);

% ==========================================
% PART 2: LOW-PASS FILTERING (3 Hz)
% ==========================================
fs = 100;       
fc = 3;         
[b, a] = butter(2, fc/(fs/2), 'low');

% Filter both the Acceleration and the Roll to remove backpack jitter
accel_mag = sqrt(imuData.AccelX.^2 + imuData.AccelY.^2 + imuData.AccelZ.^2);
accel_mag_filtered = filtfilt(b, a, accel_mag);
roll_sway_filtered = filtfilt(b, a, roll_sway);

% ==========================================
% PART 3: STEP DETECTION (IMPACTS VS SWAY)
% ==========================================
% METHOD A: Steps from Acceleration Impacts
[pks_accel, locs_accel] = findpeaks(accel_mag_filtered, 'MinPeakDistance', 40, 'MinPeakHeight', 0.5);
step_freq_accel = 1 ./ diff(t(locs_accel));

% METHOD B: Steps from Torso Sway (Roll)
[pks_roll, locs_roll_peaks] = findpeaks(roll_sway_filtered, 'MinPeakDistance', 40, 'MinPeakHeight', 1.0);
[vals_roll, locs_roll_valleys] = findpeaks(-roll_sway_filtered, 'MinPeakDistance', 40, 'MinPeakHeight', 1.0);
vals_roll = -vals_roll; 

% Combine both sides into a single timeline of steps
locs_roll_all = sort([locs_roll_peaks; locs_roll_valleys]);
step_freq_roll = 1 ./ diff(t(locs_roll_all));

% Console Output
fprintf('--- STEP DETECTION RESULTS ---\n');
fprintf('From Impacts (Accel): %d steps | Avg Freq: %.2f Hz\n', length(locs_accel), mean(step_freq_accel));
fprintf('From Sway (Roll):     %d steps | Avg Freq: %.2f Hz\n', length(locs_roll_all), mean(step_freq_roll));

% ==========================================
% PART 4: PLOTTING & LINKING AXES
% ==========================================
% Increased height to 1200 to accommodate 5 distinct subplots
figure('Name', 'BNO085 Gait Analysis', 'Color', 'w', 'Position', [100 100 900 1200]);

% Subplot 1: The Turns (Yaw)
ax1 = subplot(5,1,1);
plot(t, heading_plot, 'b', 'LineWidth', 1.5);
title('Heading Direction (Yaw)');
ylabel('Degrees');
ylim([-180 180]);
yticks([-180 -90 0 90 180]);
grid on;

% Subplot 2: Step Jolt (Magnitude)
ax2 = subplot(5,1,2);
plot(t, accel_mag_filtered, 'k', 'LineWidth', 1);
hold on;
% Solid red dots on black line
plot(t(locs_accel), pks_accel, 'r.', 'MarkerSize', 15, 'DisplayName', 'Impacts');
title('Impact Detection (Acceleration Magnitude)');
ylabel('Mag (m/s^2)');
legend('Location', 'best');
grid on;

% Subplot 3: Step Frequency from Acceleration
ax3 = subplot(5,1,3);
% Black line with dots matching the impact plot style
plot(t(locs_accel(2:end)), step_freq_accel, '-k.', 'LineWidth', 1.5, 'MarkerSize', 12);
hold on;
yline(mean(step_freq_accel), 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Avg: %.2f Hz', mean(step_freq_accel)));
title('Frequency from Impacts');
ylabel('Hz');
ylim([0 4]); 
grid on;

% Subplot 4: Torso Sway (Roll)
ax4 = subplot(5,1,4);
plot(t, roll_sway_filtered, 'Color', [0.8 0.4 0], 'LineWidth', 1.5);
hold on;
% Using red and blue dots instead of triangles
plot(t(locs_roll_peaks), pks_roll, 'r.', 'MarkerSize', 15, 'DisplayName', 'Right Shift');
plot(t(locs_roll_valleys), vals_roll, 'b.', 'MarkerSize', 15, 'DisplayName', 'Left Shift');
title('Sway Detection (Roll) with Left/Right Alternation');
ylabel('Degrees');
legend('Location', 'best');
grid on;

% Subplot 5: Step Frequency from Sway
ax5 = subplot(5,1,5);
% Orange line with dots matching the sway plot style
plot(t(locs_roll_all(2:end)), step_freq_roll, '-', 'Color', [0.8 0.4 0], 'Marker', '.', 'MarkerSize', 12, 'LineWidth', 1.5);
hold on;
yline(mean(step_freq_roll), 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Avg: %.2f Hz', mean(step_freq_roll)));
title('Frequency from Sway');
xlabel('Time (seconds)');
ylabel('Hz');
ylim([0 4]); 
grid on;

% Link the X-axes of all five subplots
linkaxes([ax1, ax2, ax3, ax4, ax5], 'x');