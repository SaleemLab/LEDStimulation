% -------------------------------------------------------------------------
% MATLAB script to simulate 60Hz PWM stimulus with 100ms raised cosine fade
% -------------------------------------------------------------------------

% 1. Define Time and Frequency parameters
Fs = 500e3;              % 500 kHz sampling rate (high enough to resolve 31.25kHz PWM)
t = 0:(1/Fs):(4 - 1/Fs); % 4 seconds total (1s off, 2s on, 1s off)

f_pwm = 31250;           % PWM frequency in Hz (matches desiredPWMFrequency)
f_sine = 60;             % Stimulus frequency in Hz
fade_duration = 0;     % 100 ms fade-in/out (matches FADE_DURATION_MS)

% 2. Generate the Stimulus Envelope (Raised Cosine Fade)
% 1s OFF (50%), 2s ON (stimulus), 1s OFF (50%)
env = zeros(size(t));
stim_start = 1;
stim_end = 3;

for i = 1:length(t)
    if t(i) >= stim_start && t(i) < (stim_start + fade_duration)
        % Fade in (Raised Cosine)
        phase = pi * (t(i) - stim_start) / fade_duration;
        env(i) = (1 - cos(phase)) / 2;
    elseif t(i) >= (stim_start + fade_duration) && t(i) <= (stim_end - fade_duration)
        % Full stimulus max amplitude
        env(i) = 1;
    elseif t(i) > (stim_end - fade_duration) && t(i) <= stim_end
        % Fade out (Raised Cosine)
        phase = pi * (t(i) - (stim_end - fade_duration)) / fade_duration;
        env(i) = (1 + cos(phase)) / 2;
    end
end

% 3. Calculate Duty Cycle 
% 0.5 (MidLumi) base, modulated by the envelope and the sine wave
duty_cycle = 0.5 + (0.5 * env .* sin(2 * pi * f_sine * t));

% 4. Generate Phase-Correct PWM (Triangle carrier)
% sawtooth(..., 0.5) generates a symmetric triangle wave (-1 to 1)
carrier = 0.5 * sawtooth(2 * pi * f_pwm * t, 0.5) + 0.5;
pwm_sig = double(duty_cycle > carrier);

% 5. Low-pass Filter to recover the 60 Hz sine wave
fc = 200; % Cutoff frequency 200 Hz (keeps 60Hz, removes 31.25kHz)
[b, a] = butter(2, fc / (Fs / 2), 'low');
pwm_filtered = filtfilt(b, a, pwm_sig);

% -------------------------------------------------------------------------
% Plotting
% -------------------------------------------------------------------------
figure('Name', '60Hz Sine PWM Stimulus Simulation', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Top Plot: Raw PWM Signal
ax1 = subplot(2, 1, 1);
plot(t, pwm_sig, 'b', 'LineWidth', 0.00001);
title('Raw PWM Signal (31.25 kHz Carrier)');
xlabel('Time (s)');
ylabel('Logic Level (0 / 1)');
ylim([-0.1 1.1]);

% Bottom Plot: Low-pass Filtered Signal
ax2 = subplot(2, 1, 2);
plot(t, pwm_filtered, 'r', 'LineWidth', 1.5);
title('Low-Pass Filtered Signal (Recovered 60 Hz Sine Wave)');
xlabel('Time (s)');
ylabel('Analog Equivalent (Duty Cycle)');
ylim([-0.1 1.1]);
grid on;

% Link the X-axes so zooming on one zooms on the other
linkaxes([ax1, ax2], 'x');

disp('Simulation complete! Note: Zoom in heavily on the X-axis of the top plot to see individual PWM pulses.');