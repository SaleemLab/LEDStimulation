function [totalDurationSec, freqList, timeList] = calculateSteppedSweepDuration(startFreq, endFreq, stepFreq, cyclesPerFreq)
    % ---------------------------------------------------------------------
    % Calculates the total duration of the SteppedFrequencySweep stimulus.
    % 
    % INPUTS:
    %   startFreq     : Starting frequency in Hz
    %   endFreq       : Ending frequency in Hz
    %   stepFreq      : Amount to step the frequency by (automatically made positive)
    %   cyclesPerFreq : Number of complete sine waves per frequency step
    %
    % OUTPUTS:
    %   totalDurationSec : Total time in seconds
    %   freqList         : Array of the exact frequencies hit during the sweep
    %   timeList         : Array of how much time was spent at each frequency
    % ---------------------------------------------------------------------

    % Match C++ logic: stepFreq is forced positive
    stepFreq = abs(stepFreq);
    sweepingUp = startFreq <= endFreq;
    
    currentFreq = startFreq;
    totalDurationSec = 0;
    
    % Arrays to hold the data for plotting
    freqList = [];
    timeList = [];
    
    % Exact replica of the C++ while loop logic
    while (sweepingUp && currentFreq <= endFreq) || (~sweepingUp && currentFreq >= endFreq)
        
        % Safety check to prevent divide-by-zero
        if currentFreq <= 0
            warning('Frequency reached 0. Aborting further calculation to prevent Infinity.');
            break;
        end
        
        % Calculate time spent at this specific frequency step
        timeAtFreq = cyclesPerFreq / currentFreq;
        totalDurationSec = totalDurationSec + timeAtFreq;
        
        % Store for tracking
        freqList(end+1) = currentFreq;
        timeList(end+1) = timeAtFreq;
        
        % Step the frequency
        if sweepingUp
            currentFreq = currentFreq + stepFreq;
        else
            currentFreq = currentFreq - stepFreq;
        end
    end
    
    % Account for the hardware delay at the end of the Arduino function
    % delayMicroseconds(2000); -> 0.002 seconds
    totalDurationSec = totalDurationSec + 0.002;
    
    % --- Command Window Output ---
    fprintf('===================================================\n');
    fprintf('Stepped Sweep Parameters:\n');
    fprintf('Freq Range : %.2f Hz to %.2f Hz (Step: %.2f Hz)\n', startFreq, endFreq, stepFreq);
    fprintf('Cycles/Step: %d\n', cyclesPerFreq);
    fprintf('Total Steps: %d\n', length(freqList));
    fprintf('---------------------------------------------------\n');
    fprintf('TOTAL DURATION: %.4f seconds (%.1f ms)\n', totalDurationSec, totalDurationSec * 1000);
    fprintf('===================================================\n\n');
    
    % --- Optional: Plot the Time Distribution ---
    figure('Name', 'Stepped Sweep Duration Analysis', 'Position', [100, 100, 800, 400]);
    bar(freqList, timeList, 'FaceColor', [0.2 0.6 0.8]);
    xlabel('Frequency (Hz)');
    ylabel('Time Spent at Frequency (Seconds)');
    title(sprintf('Time Distribution per Frequency Step (Total: %.2fs)', totalDurationSec));
    grid on;
end