%% Pupil Labs Neon IMU: Orientation & Step Rate Analysis
clear; clc; close all;

% 1. Load Data
imu = readtable('imu.csv');

% 2. Clean Variable Names
imu.Properties.VariableNames = {...
    'SectionID', 'RecordingID', 'Time_ns', ...
    'GX', 'GY', 'GZ', ...
    'AX', 'AY', 'AZ', ...
    'Roll', 'Pitch', 'Yaw', ...
    'QW', 'QX', 'QY', 'QZ'};

% 3. Calculate Time and Vertical Filter
time_s = (imu.Time_ns - imu.Time_ns(1)) / 1e9;
fs = 1 / mean(diff(time_s)); % Sampling Frequency

% Isolate Vertical Axis (assuming AZ is gravity-aligned)
% We subtract the mean (1g) to center the "bounce" around 0
v_accel = imu.AZ - mean(imu.AZ); 

% Low-pass filter at 3Hz to isolate the walking rhythm
[b, a] = butter(4, 3/(fs/2), 'low');
clean_v = filtfilt(b, a, v_accel);

% 4. Frequency Analysis (FFT) for Step Rate
L = length(clean_v);
Y = fft(clean_v);
P2 = abs(Y/L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

% Find the dominant frequency in the human walking range (0.5 - 4 Hz)
walk_range = (f > 0.5 & f < 4);
[peak_amp, peak_idx] = max(P1(walk_range));
f_walk = f(walk_range);
step_freq_hz = f_walk(peak_idx);

% 5. Integrated Visualization
figure('Color', 'w', 'Position', [100, 100, 900, 900]);
t = tiledlayout(2,1, 'TileSpacing', 'compact');
title(t, ['Walk-Turn-Walk Analysis | Detected Step Rate: ' num2str(step_freq_hz, '%.2f') ' Hz']);

% --- Plot 1: Orientation (The Path) ---
nexttile
plot(time_s, imu.Yaw, 'LineWidth', 2, 'Color', '#D95319')
ylabel('Yaw (deg)')
title('Heading (Orientation)')
grid on
legend('Yaw / Heading')

% --- Plot 2: The "Bounce" (The Steps) ---
nexttile
plot(time_s, v_accel, 'Color', [0.8 0.8 0.8]) % Raw in gray
hold on
plot(time_s, clean_v, 'b', 'LineWidth', 1.5) % Filtered in blue
ylabel('Vertical Accel (g)')
title('Gait Rhythm (Vertical Oscillation)')
grid on
legend('Raw AZ', 'Filtered Steps')

xlabel(t, 'Time (seconds)')