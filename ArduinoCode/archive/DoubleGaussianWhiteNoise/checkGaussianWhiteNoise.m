% --- Script to Test for Gaussian White Noise (No Econometrics Toolbox - Refined ACF Check - Combined Plots) ---

clear;
clc;
close all; % Close any existing figures

% --- Prerequisites ---
% Requires: Statistics and Machine Learning Toolbox 
% Optional: Signal Processing Toolbox 
fprintf('Checking required/optional toolboxes...\n');
if ~license('test', 'Statistics_Toolbox')
    error('This script requires the Statistics and Machine Learning Toolbox.');
else
    fprintf('- Statistics and Machine Learning Toolbox: Found.\n');
end
hasSignalToolbox = license('test', 'Signal_Toolbox');
if hasSignalToolbox
    fprintf('- Signal Processing Toolbox: Found (will use its autocorr).\n');
else
    fprintf('- Signal Processing Toolbox: Not found (will calculate ACF manually).\n');
end

% --- 1. Generate Sample Data ---
rng('shuffle'); % Use a different seed each time
fprintf('\nUsing rng(''shuffle'') for potentially different random sequence.\n');

% N = 1000; % Number of data points
% mu = 0; sigma = 2;
% x_gw = mu + sigma * randn(N, 1); 
% x = x_gw; % Test the generated Gaussian white noise
% fprintf('Testing data with N = %d samples.\n', N);

x = [-2.29,0.28,-0.41,-1.62,-1.62,2.50,0.97,0.65,-0.45,-0.48,1.96,-0.49,0.39,-0.27,1.31,-1.40,1.08,0.59,-0.43,0.77,0.34,-0.46,-0.05,-0.20,-0.19,0.99,-1.09,-0.17,-0.21,0.57,0.21,0.19,-0.94,-1.04,1.01,-1.10,-1.09,-0.43,0.92,-0.65,-1.25,0.30,-0.74,1.88,-0.16,-1.92,0.04,-0.20,-0.26,0.72,0.39,-0.44,-0.55,1.02,0.42,-0.74,1.48,1.43,0.34,0.57,-0.23,-1.16,-0.82,0.65,-0.10,-0.95,-0.21,0.62,1.36,-1.03,1.89,0.47,0.94,0.74,0.40,-0.45,0.07,-0.13,0.61,1.23,-0.28,-0.85,0.02,1.99,-1.69,0.47,-0.39,0.47,0.48,0.56,-0.28,-1.67,-1.60,0.04,0.53,0.15,-0.05,1.28,1.54,-1.08,0.29,1.05,-0.04,1.17,-1.00,-0.38,0.37,-0.09,1.05,-0.60,1.02,0.89,-0.69,-0.03,-1.46,-0.62,-0.83,0.35,-0.07,-0.43,-1.14,-1.04,-0.92,0.20,1.12,0.45,0.47,0.12,1.18,-1.78,-0.38,0.80,-0.34,-0.23,-1.75,0.87,-2.01,-0.61,-1.52,-0.31,-0.10,-1.16,-1.12,0.05,0.85,1.10,-0.92,-0.37,0.31,0.88,1.39,1.49,1.45,0.77,-1.22,-0.77,0.60,-1.41,-0.94,1.10,-2.54,-0.18,0.36,0.65,0.21,-0.53,-0.87,0.76,-1.39,0.18,-0.42,-0.00,-1.04,1.35,-0.05,-0.27,0.62,1.23,-0.27,-0.46,0.14,0.47,-0.46,-1.68,-1.23,0.44,-0.86,-0.98,0.00,0.35,-0.32,-1.15,0.16,-0.17,0.16,0.61,0.41,-0.34,0.28,-1.03,0.16,-1.19,0.81,0.55,0.62,-0.96,-0.80,0.62,-0.55,0.39,-0.41,0.05,1.26,0.21,0.09,-0.38,-1.37,-0.38,-1.84,0.19,0.93,-0.82,-0.62,0.49,1.65,1.37,-0.31,0.44,-0.42,-0.67,0.26,0.10,-1.47,-0.34,0.32,-2.45,-0.36,-0.57,1.36,-0.24,1.37,-1.44,-0.15,0.24,1.04,0.93,1.18,-0.23,-0.50,-0.10,-0.08,0.05,-0.03,0.92,0.46,-0.15,-0.81,0.27,-0.14,0.09,0.19,0.58,-0.26,-1.58,-1.13,1.29,-1.21,-0.46,-0.77,-0.52,-0.15,0.99,0.27,-0.08,-0.93,-2.34,-0.60,-0.53,-2.17,0.33,-1.10,0.80,0.58,-0.65,0.30,-0.08,0.44,-0.65,0.82,-3.42,-0.54,-1.19,2.76,0.74,-1.57,-0.15,-0.20,0.33,-0.61,0.68,-1.09,-0.06,0.46,-0.16,-0.14,-0.96,1.85,1.42,0.38,-0.32,-0.00,-0.17,-1.50,0.56,0.89,2.55,-0.06,-1.60,0.13,-1.23,0.41,-1.41,0.58,1.55,-0.08,-0.13,-2.25,1.19,0.93,0.39,-0.82,-0.06,-0.35,0.27,-0.92,0.54,1.35,-0.26,1.55,1.88,0.99,-1.13,2.14,1.42,0.24,-1.42,0.10,-0.28,1.47,-0.36,-0.15,0.30,0.36,-0.14,0.13,0.12,0.09,-0.53,-0.55,1.82,-0.42,-0.04,-0.56,-0.87,0.29,0.09,0.21,0.81,-0.76,0.10,-0.71,1.92,-0.15,1.71,0.27,-0.39,1.13,0.26,0.33,0.81,-0.22,0.15,0.72,-1.15,0.67,-0.38,-1.27,0.22,2.04,-1.59,0.96,0.58,-1.11,-0.05,1.73,-0.50,1.32,-0.19,0.49,0.68,1.00,0.28,-1.25,1.46,-0.99,-0.57,-0.29,-1.11,2.94,1.47,-0.16,-0.39,0.60,-1.84,-1.48,-0.75,1.68,-0.31,-0.14,-1.33,-0.30,-0.14,-0.59,-0.18,-2.13,0.01,-0.53,-0.22,-0.57,0.46,-1.51,-0.49,1.43,-0.19,-0.28,-2.88,0.53,0.04,0.17,-2.73,-1.22,-1.88,-0.57,-0.04,-1.01,1.20,1.24,0.02,0.77,-0.65,-1.65,-1.19,-0.34,-0.71,-0.11,-1.14,-1.49,1.17,1.04,3.09,-1.11,-2.31,-2.28,1.09,-0.46,-0.04,0.38,-1.30,-0.40,0.35,-0.10,-2.38,-1.28,-0.65,-0.28,0.67,-0.62,0.92,-1.23,0.11,-1.26,1.02,2.03,0.84,0.79,-0.52,-0.70,-0.88,-0.19,-0.54,0.45,2.49,-0.83,0.96,-1.00,2.37,0.16,-0.12,0.32,0.72,0.95,0.94,-0.42,0.44,0.55,-0.24,-0.42,-0.13,0.16,0.26,0.07,0.33,-0.59,-0.42,-1.12,-1.46,0.56,-0.24,-1.05,0.85,-0.09,0.57,1.06,-0.66,-1.21,-2.50,2.32,-1.48,1.41,-1.00,-0.04,0.97,1.44,-1.06,-0.86,-2.01,1.14,3.05,0.50,-0.27,0.07,-0.80,0.22,-0.20,-0.76,-0.21,-0.20,0.92,1.52,-0.62,1.40,-0.24,-1.96,-2.36,0.24,-1.64,-0.78,0.91,-0.84,-2.19,-0.85,1.37,0.80,-2.26,2.46,0.40,-2.07,-1.88,-0.05,-0.55,-1.29,-1.19,-0.52,0.27,-1.53,0.98,2.73,0.10,0.79,1.64,-0.72,-1.24,-0.83,-2.73,0.19,0.94,-1.70,-1.47,1.81,2.23,1.02,-0.53,-1.09,-1.17,2.40,-1.07,-0.03,0.27,0.10,-1.26,0.29,0.19,-0.94,0.68,-0.03,2.02,-0.53,1.11,-0.85,0.29,0.37,1.14,0.96,-0.28,-1.10,-0.48,-1.74,0.22,-0.93,-0.76,-0.37,1.40,-1.80,-0.24,1.15,2.00,-0.53,-0.85,1.39,1.60,-1.24,-1.11,2.60,-0.75,1.61,1.98,0.03,1.29,0.15,0.91,0.23,-0.29,0.39,-0.40,0.24,-0.47,1.19,-0.13,0.20,-0.20,-0.26,-0.50,-0.65,-0.25,1.48,1.11,0.71,1.14,1.83,-0.05,-0.63,-0.36,-1.64,-0.82,0.09,-1.16,1.67,0.50,2.33,0.19,-0.67,0.01,-1.38,-0.77,-0.76,-0.02,-1.48,0.57,-0.08,0.19,-1.81,-2.05,-1.29,0.23,0.63,0.43,-0.36,-1.64,-1.13,-1.20,-0.11,0.21,0.66,-1.13,-0.33,-0.40,0.90,2.02,-1.03,-1.34,-0.40,-0.59,1.91,0.78,-0.44,-0.65,0.07,1.37,-0.89,0.24,0.66,0.30,-1.33,-0.07,-1.56,-0.62,-0.43,-0.57,-0.16,-1.05,0.52,-1.45,0.47,0.61,1.13,-1.28,0.52,0.70,0.91,0.48,0.92,-0.29,-0.16,-0.02,-0.51,0.51,-0.48,-1.03,0.30,-0.02,-0.47,0.44,0.65,-0.69,0.20,1.38,-0.80,-0.25,-0.93,-0.28,-0.19,-0.10,-2.23,2.52,-0.86,2.50,1.92,0.28,0.54,0.51,0.54,-0.11,-0.88,-0.05,-0.22,-0.83,-1.32,-1.46,1.12,-0.48,0.08,0.09,1.63,-0.57,0.50,-2.80,1.75,0.38,-1.06,1.65,0.10,0.47,-0.89,0.63,-0.62,0.50,0.74,1.55,0.29,0.42,0.68,-0.43,0.26,0.47,0.49,-0.41,0.29,0.52,-1.09,-1.52,-0.24,-1.98,-0.51,1.04,-1.03,-1.14,0.77,1.84,-0.19,1.78,1.01,1.56,1.36,-0.52,-0.71,-0.57,0.24,-0.88,-1.26,1.76,0.77,-1.35,-0.16,0.23,0.40,-0.13,-0.11,0.34,-0.50,1.07,-1.72,-1.88,-1.12,-0.12,-0.14,0.95,-1.27,-0.49,1.26,1.56,-0.03,0.91,-2.31,-0.85,-0.35,1.75,-0.68,-0.74,-0.23,0.94,1.00,1.97,0.84,0.82,-0.24,0.17,0.57,1.32,1.00,-1.32,0.08,0.72,1.16,-0.36,0.97,0.14];
N = numel(x);

% --- 2. Set Significance Level ---
alpha = 0.05; 
fprintf('Significance level alpha = %.2f\n', alpha);

% --- 3. Perform Statistical Tests ---

% --- Test 3.1: Zero Mean (Requires Statistics Toolbox) ---
fprintf('\n--- Test 1: Zero Mean (t-test) ---\n');
[h_mean, p_mean, ci_mean, ~] = ttest(x, 0, 'Alpha', alpha); 
mean_ok = (h_mean == 0);
fprintf(' H0: Mean is zero. Result: %s H0. Mean is %s with zero (p=%.4f).\n', ...
    iif(mean_ok,'Fail to reject','Reject'), iif(mean_ok,'statistically consistent','significantly different from'), p_mean);

% --- Test 3.2: Normality (Requires Statistics Toolbox) ---
fprintf('\n--- Test 2a: Normality (Jarque-Bera) ---\n');
[h_jb, p_jb, ~, ~] = jbtest(x, alpha);
normality_jb_ok = (h_jb == 0);
fprintf(' H0: Data is normally distributed. Result: %s H0 (p=%.4f).\n', ...
    iif(normality_jb_ok,'Fail to reject','Reject'), p_jb);

fprintf('\n--- Test 2b: Normality (Lilliefors) ---\n');
[h_lillie, p_lillie, ~, ~] = lillietest(x, 'Alpha', alpha);
normality_lillie_ok = (h_lillie == 0);
fprintf(' H0: Data is normally distributed. Result: %s H0 (p=%.4f).\n', ...
     iif(normality_lillie_ok,'Fail to reject','Reject'), p_lillie);
% Decide if both or just one needs to pass for overall normality check later
normality_ok = normality_jb_ok; % Using JB test result primarily


% --- Test 3.3: Uncorrelated (Independence) ---
fprintf('\n--- Test 3: Uncorrelated (Autocorrelation & Randomness Checks) ---\n');
maxLag = min(30, floor(N/4)); 
fprintf(' Checking ACF up to lag %d.\n', maxLag);



    x_centered = x - mean(x); 
    [acf_raw, lags_raw] = xcorr(x_centered, maxLag, 'coeff'); 
    zeroLagIndex = find(lags_raw == 0);
    acf_vals = acf_raw(zeroLagIndex:end);
    lags = lags_raw(zeroLagIndex:end);
    bound_val = 1.96 / sqrt(N); 
    bounds = [bound_val; -bound_val]; 
    fprintf(' Approx. 95%% confidence bound magnitude: +/- %.4f\n', bound_val);


% --- Refined Significance Check for ACF ---
lags_to_check = lags(2:end);          % Lags > 0
acf_vals_to_check = acf_vals(2:end);  % Corresponding ACF values
significant_indices = find(abs(acf_vals_to_check) > bound_val);
significant_lags = lags_to_check(significant_indices);
num_significant = length(significant_lags);
allowed_significant = max(1, floor(alpha * maxLag)); 

fprintf(' ACF Check:\n');
fprintf(' - Found %d significant lag(s) out of %d checked (lags 1-%d).\n', ...
        num_significant, maxLag, maxLag);
if num_significant > 0
    fprintf(' - Significant lag(s): %s\n', num2str(significant_lags'));
end
fprintf(' - Tolerance for statistical chance (approx %.0f%% of %d): %d lag(s).\n', ...
        alpha*100, maxLag, allowed_significant);
autocorr_ok = (num_significant <= allowed_significant); 
if autocorr_ok
    fprintf(' - Result: ACF check PASSED (number of significant lags within tolerance).\n');
else
    fprintf(' - Result: ACF check FAILED (number of significant lags exceeds tolerance).\n');
    fprintf('   (Visual inspection recommended to assess the pattern of correlation.)\n');
end

% Formal Test for Randomness: Runs Test (Requires Statistics Toolbox)
fprintf('\n Runs Test (Test for independence/randomness):\n');
[h_run, p_run, ~] = runstest(x, alpha); 
runs_ok = (h_run == 0);
fprintf(' H0: Sequence is random. Result: %s H0 (p=%.4f).\n', ...
    iif(runs_ok,'Fail to reject','Reject'), p_run);

% Overall assessment for correlation/independence:
uncorrelated_ok = autocorr_ok && runs_ok; 
fprintf('\n Overall Independence Result: %s\n', iif(uncorrelated_ok, 'PASSES', 'FAILS'));


% --- 4. Visualizations (Combined Figure) ---
fprintf('\nGenerating combined plot figure...\n');
figure; % Create the single figure window

% Time series plot
subplot(2, 2, 1); 
plot(x); 
title('Time Series Plot'); xlabel('Sample Index'); ylabel('Value'); 
grid on; axis tight;

% Histogram with normal fit
subplot(2, 2, 2); 
histfit(x); 
title('Histogram with Normal Fit'); xlabel('Value'); ylabel('Frequency'); 
grid on;

% Q-Q plot 
subplot(2, 2, 3); 
qqplot(x); 
title('Q-Q Plot vs Standard Normal'); 
grid on;

% ACF Plot (Moved to subplot 2,2,4)
subplot(2, 2, 4); 
stem(lags, acf_vals, 'filled', 'MarkerSize', 4);
hold on;
plot(lags, bounds(1)*ones(size(lags)), 'r--', 'LineWidth', 1.5); 
plot(lags, bounds(2)*ones(size(lags)), 'r--', 'LineWidth', 1.5); 
plot(lags, zeros(size(lags)), 'k-'); 
% Highlight significant lags found
if num_significant > 0
   plot(significant_lags, acf_vals(significant_indices+1), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); 
   legend('ACF', 'Conf. Bounds', '', 'Sig. Lags', 'Location', 'northeast'); % Adjusted legend location potentially
else
   legend('ACF', 'Conf. Bounds', 'Location', 'northeast'); 
end
hold off;
xlabel('Lag'); ylabel('ACF');
title(sprintf('Sample ACF (Lags 0-%d)', maxLag));
grid on; ylim([-1 1]); xlim([0 maxLag]);

% Adjust overall figure appearance if needed
sgtitle('Gaussian White Noise Test Diagnostics'); % Optional super title


% --- 5. Overall Conclusion ---
fprintf('\n--- Overall Conclusion (alpha=%.2f) ---\n', alpha);

is_gaussian_white_noise = mean_ok && normality_ok && uncorrelated_ok; 

if is_gaussian_white_noise
    fprintf('Verdict: The data sequence IS statistically consistent with Gaussian White Noise.\n');
    fprintf('(Based on t-test, normality test(s), tolerant ACF check, and Runs test).\n');
else
    fprintf('Verdict: The data sequence IS NOT statistically consistent with Gaussian White Noise.\n');
    fprintf('(Based on t-test, normality test(s), tolerant ACF check, and Runs test).\n');
    fprintf('Reason(s):\n');
    if ~mean_ok
        fprintf(' - The mean test failed.\n');
    end
    if ~normality_ok 
        fprintf(' - The normality test(s) failed.\n');
    end
    if ~uncorrelated_ok
        fprintf(' - The independence checks failed:\n');
         if ~autocorr_ok
            fprintf('   - Too many significant lags found in ACF (> %d).\n', allowed_significant);
         end
         if ~runs_ok
             fprintf('   - The sequence failed the Runs test for randomness.\n');
         end
    end
end

% Helper function
function out = iif(condition, true_val, false_val)
    if condition, out = true_val; else, out = false_val; end
end