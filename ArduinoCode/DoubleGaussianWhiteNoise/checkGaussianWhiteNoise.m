% --- Script to Test for Gaussian White Noise (No Econometrics Toolbox - Refined ACF Check - Combined Plots) ---

% clear;
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

% x = vals;
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