% --- Parameters Definition ---

fs = 1000;                  % Sampling frequency in Hz
totalDuration = 10;         % Total duration of the plot in seconds
stimDuration = 2;           % Duration of each stimulus in seconds
isiDuration = 2;            % Duration of the inter-stimulus interval (mean luminance) in seconds
meanLuminance = 0.5;        % Mean luminance level (e.g., on a 0 to 1 scale)
amplitude = 0.4;            % Amplitude of the sinewave, must be <= meanLuminance to stay within [0,1]

% --- Stimulus properties ---
flickerFrequencies = [2, 10, 40]; % Possible frequencies for the sinewave in Hz

% --- Time Vector ---
t = 0:1/fs:totalDuration-1/fs; % Time vector for the entire duration

% --- Initialize Signal Vectors ---
% Start with both signals at mean luminance
greenSignal = ones(size(t)) * meanLuminance;
magentaSignal = ones(size(t)) * meanLuminance;

% --- Generate the Stimulus Sequence using a while loop ---
currentTime = 0; % Keeps track of the start time for the next stimulus block

fprintf('--- Stimulus Generation Log ---\n');

% Loop as long as the start of the next stimulus is within the total duration
while currentTime < totalDuration
    
    % --- Define start and end times for the current stimulus ---
    stimStartTime = currentTime;
    stimEndTime = stimStartTime + stimDuration;

    % --- Randomly select properties for the current stimulus ---
    currentFreq = flickerFrequencies(randi(length(flickerFrequencies)));
    flickeringColor = randi([1, 2]); % 1 for Green, 2 for Magenta

    % --- Define start and end points in array indices ---
    % Ensure the stimulus does not exceed the total duration
    startIdx = round(stimStartTime * fs) + 1;
    endIdx = round(min(stimEndTime, totalDuration) * fs);

    % --- Create the sinewave for the current stimulus ---
    currentStimDuration = (endIdx - startIdx + 1) / fs;
    stimulusTime = 0:1/fs:currentStimDuration-1/fs;
    sinewave = amplitude * sin(2 * pi * currentFreq * stimulusTime) + meanLuminance;

    % --- Assign the sinewave to the selected color channel ---
    if flickeringColor == 1 % Green flickers
        greenSignal(startIdx:endIdx) = sinewave;
        fprintf('Time %.1f-%.1fs: Green stimulus at %d Hz\n', stimStartTime, stimEndTime, currentFreq);
    else % Magenta flickers
        magentaSignal(startIdx:endIdx) = sinewave;
        fprintf('Time %.1f-%.1fs: Magenta stimulus at %d Hz\n', stimStartTime, stimEndTime, currentFreq);
    end

    % --- Update currentTime for the start of the next cycle ---
    currentTime = currentTime + stimDuration + isiDuration;
end

% --- Plotting the Stimuli ---
figure; % Create a new figure window
hold on; % Hold on to plot multiple lines on the same axes

plot(t, greenSignal, 'g', 'LineWidth', 2);
plot(t, magentaSignal, 'm', 'LineWidth', 2);

hold off; % Release the plot hold

% --- Formatting the Plot ---
title('Example Sinewave Stimuli (10s Segment)');
xlabel('Time (s)');
ylabel('Luminance');
ylim([-0.1, 1.1]); % Set y-axis limits to view the full wave
legend('Green Channel', 'Magenta Channel', 'Location', 'northeast');
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 12);

%%

% --- Parameters Definition ---

fs = 1000;                  % Sampling frequency in Hz
totalDuration = 10;         % Total duration of the plot in seconds
stimDuration = 2;           % Duration of each stimulus in seconds
isiDuration = 2;            % Duration of the inter-stimulus interval (mean luminance) in seconds
meanLuminance = 0.5;        % Mean luminance level (e.g., on a 0 to 1 scale)

% --- Stimulus Properties ---
flickerFrequencies = [2, 10, 40]; % Possible frequencies in Hz
contrastLevels = [0.05, 0.10, 0.20, 0.40]; % Possible contrasts (5%, 10%, 20%, 40%)

% --- Time Vector ---
t = 0:1/fs:totalDuration-1/fs; % Time vector for the entire duration

% --- Initialize Signal Vectors ---
% Start with both signals at mean luminance
greenSignal = ones(size(t)) * meanLuminance;
magentaSignal = ones(size(t)) * meanLuminance;

% --- Generate the Stimulus Sequence using a while loop ---
currentTime = 0; % Keeps track of the start time for the next stimulus block

fprintf('--- Stimulus Generation Log ---\n');

% Loop as long as the start of the next stimulus is within the total duration
while currentTime < totalDuration
    
    % --- Define start and end times for the current stimulus ---
    stimStartTime = currentTime;
    stimEndTime = stimStartTime + stimDuration;

    % --- Randomly select properties for the current stimulus ---
    currentFreq = flickerFrequencies(randi(length(flickerFrequencies)));
    currentContrast = contrastLevels(randi(length(contrastLevels)));
    flickeringColor = randi([1, 2]); % 1 for Green, 2 for Magenta

    % --- Calculate the amplitude based on the selected contrast ---
    amplitude = currentContrast * meanLuminance;

    % --- Define start and end points in array indices ---
    startIdx = round(stimStartTime * fs) + 1;
    endIdx = round(min(stimEndTime, totalDuration) * fs);

    % --- Create the sinewave for the current stimulus ---
    currentStimDuration = (endIdx - startIdx + 1) / fs;
    stimulusTime = 0:1/fs:currentStimDuration-1/fs;
    sinewave = amplitude * sin(2 * pi * currentFreq * stimulusTime) + meanLuminance;

    % --- Assign the sinewave to the selected color channel ---
    if flickeringColor == 1 % Green flickers
        greenSignal(startIdx:endIdx) = sinewave;
        fprintf('Time %.1f-%.1fs: Green stimulus at %d Hz, %.0f%% contrast\n', ...
            stimStartTime, stimEndTime, currentFreq, currentContrast*100);
    else % Magenta flickers
        magentaSignal(startIdx:endIdx) = sinewave;
        fprintf('Time %.1f-%.1fs: Magenta stimulus at %d Hz, %.0f%% contrast\n', ...
            stimStartTime, stimEndTime, currentFreq, currentContrast*100);
    end

    % --- Update currentTime for the start of the next cycle ---
    currentTime = currentTime + stimDuration + isiDuration;
end

% --- Plotting the Stimuli ---
figure; % Create a new figure window
hold on; % Hold on to plot multiple lines on the same axes

plot(t, greenSignal, 'g', 'LineWidth', 1);
plot(t, magentaSignal, 'm', 'LineWidth', 1);

hold off; % Release the plot hold

% --- Formatting the Plot ---
title('Example Sinewave Stimuli with Variable Contrast');
xlabel('Time (s)');
ylabel('Luminance');
ylim([-0.1, 1.1]); % Set y-axis limits to view the full wave
legend('Green Channel', 'Magenta Channel', 'Location', 'northeast');
grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 12);