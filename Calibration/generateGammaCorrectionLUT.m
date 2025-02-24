%% Gamma calibration

outputFileName = 'GREEN_120mA_240225_gammaLUT.txt';

%% parse power values from gamma data

power_table = green120mA;

figure, 
subplot(121)
plot(power_table.Time,power_table.Power)

waitTime = 10000; % duration of each power level used in gamma correction procedure
delayTime = 2000; % how long to wait before taking values
takeValsDuration = 6000;
dutyCyclesOrig = 0:1:100;

%%
startTime = 11289; % estimate from figure


for icycle = 1:numel(dutyCyclesOrig)
    temp_startTime = startTime + waitTime*(icycle-1) + waitTime;
    temp_endTime = temp_startTime + takeValsDuration;
    temp_idx = find(power_table.Time>=temp_startTime,1,'first'):find(power_table.Time<=temp_endTime,1,'last');
    temp_powVals = power_table.Power(temp_idx);
    powerValsOrig(icycle) = median(temp_powVals);
end

subplot(122)
plot(dutyCyclesOrig, powerValsOrig)


%%


% Define input parameters
num_levels = 1041; % Number of possible pixel values
pixel_values = linspace(0, num_levels - 1, num_levels); % 0 to 1040

% Define your measured brightness values (Replace with actual data)
duty_cycle = dutyCyclesOrig;
measured_brightness = powerValsOrig / max(powerValsOrig); % Normalize to [0,1]

% Ensure brightness is monotonically increasing using cumulative max
monotonic_brightness = cummax(measured_brightness);

% Remove duplicate brightness values
[unique_brightness, unique_idx] = unique(monotonic_brightness, 'stable');
unique_duty_cycle = duty_cycle(unique_idx); % Keep corresponding duty cycles

% Create a target linear brightness curve
target_brightness = linspace(0, 1, num_levels); % Ideal linear brightness

% Interpolate to find corrected pixel values
corrected_pixel_values = interp1(unique_brightness, unique_duty_cycle, target_brightness, 'linear', 'extrap');

% Scale back to pixel range
LUT = round(corrected_pixel_values * (num_levels - 1) / 100);
LUT = max(0, min(num_levels - 1, LUT)); % Clip values to valid range

% Ensure LUT is strictly increasing
LUT = cummax(LUT); 

% Enforce minimum value to off
LUT(1)=0;

% Display LUT graph
figure;
plot(pixel_values, LUT, 'b', 'LineWidth', 2);
xlabel('Input Pixel Value');
ylabel('Corrected Pixel Value');
title('Monotonic Gamma Correction LUT');
grid on;



%% Load or generate your LUT
LUT = round(LUT);  % Ensure integer values
num_values = length(LUT);

% Format as a C-style array (10 values per line for readability)
fileID = fopen(outputFileName, 'w');
fprintf(fileID, 'const uint16_t LUT[%d] = {\n    ', num_values);

for i = 1:num_values
    fprintf(fileID, '%d', LUT(i));
    if i < num_values
        fprintf(fileID, ', ');  % Add comma between values
    end
    if mod(i, 10) == 0  % Break line every 10 values for readability
        fprintf(fileID, '\n    ');
    end
end

fprintf(fileID, '};\n');
fclose(fileID);

disp(['LUT saved to ', outputFileName]);