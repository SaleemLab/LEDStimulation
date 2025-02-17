function [wavelengths, power_per_bin, P_true, correction_factor] = getLEDSpectraFromPowerMeter(...
    P_total, lambda_meas, LED_wavelengths, LED_spectrum, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity)

% Raw LED spectrum (in arbitrary units)
S = LED_spectrum;

% Interpolate the power meter responsivity to the LED wavelengths
R_vals = interp1(P_meter_cal_wavelengths, P_meter_cal_responsivity, ...
    LED_wavelengths, 'spline');
% For clarity, set R_wavelengths equal to LED_wavelengths
R_wavelengths = LED_wavelengths;

% Find the responsivity at the measurement (calibration) wavelength using nearest neighbor:
[~, idx] = min(abs(R_wavelengths - lambda_meas));
R_peak = R_vals(idx);

% Compute the correction factor:
% Measured power (P_total) is given by: 
%   P_total = k/R_peak * trapz(S .* R_vals)
% True power should be: 
%   P_true = k * trapz(S)
% So, correction_factor = P_true / P_total = (R_peak * trapz(S)) / trapz(S .* R_vals)
correction_factor = R_peak * trapz(LED_wavelengths, S) / trapz(LED_wavelengths, S .* R_vals);

% Estimate the true total LED power
P_true = P_total * correction_factor;

% Now distribute P_true according to the LED spectral shape
P_lambda = P_true * S / trapz(LED_wavelengths, S);  % Units: Watts per nm

% Assuming uniform spacing, calculate power per bin (in Watts)
delta_lambda = LED_wavelengths(2) - LED_wavelengths(1);
power_per_bin = P_lambda * delta_lambda;  

% Return wavelengths for convenience
wavelengths = LED_wavelengths;

%% Plot the Results


% figure;subplot(211)
% hold on,
% yyaxis left
% plot(LED_wavelengths, S, 'b-', 'LineWidth',1.5);
% xlabel('Wavelength (nm)');
% ylabel('LED Spectrum (a.u.)');
% grid on;
% yyaxis right
% plot(LED_wavelengths, power_per_bin, 'r-', 'LineWidth',1.5);
% 
% subplot(212)
% plot(P_meter_cal_wavelengths, P_meter_cal_responsivity, 'k:')
% ylabel('Power meter responsivity')
% title(fprintf('Correction factor = %.3f\n', correction_factor))
% % Display the measured and true power values
% fprintf('Measured power (biased) = %.3f W\n', P_total);
% fprintf('Estimated true power = %.3f W\n', P_true);
% fprintf('Correction factor = %.3f\n', correction_factor);
