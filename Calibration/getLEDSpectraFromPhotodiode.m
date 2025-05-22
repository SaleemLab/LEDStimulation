function [wavelengths, power_per_bin, P_true] = getLEDSpectraFromPhotodiode(...
    V_measured, V_dark, LED_wavelengths, LED_spectrum, ...
    PD_cal_wavelengths, PD_cal_responsivity, PD_gain_value)


% Calculate voltage signal
V_signal = V_measured - V_dark;

% Calculate photocurrent (Ip)
Ip_A = V_signal / PD_gain_value; % amps

% Interpolate Photodiode Responsivity to match LED wavelength points
% This ensures R(lambda_i) is known for each S_LED_norm(lambda_i)
R_interp_AW = interp1(PD_cal_wavelengths, PD_cal_responsivity, LED_wavelengths, 'linear', 'extrap');
% Handle potential NaN from extrapolation if LED spectrum is outside photodiode range
R_interp_AW(isnan(R_interp_AW)) = 0; % Assume zero responsivity outside defined range

% Calculate Effective Responsivity (R_eff)
% Since bin size is 1 nm, delta_lambda terms cancel out.
numerator_Reff = sum(LED_spectrum .* R_interp_AW);
denominator_Reff = sum(LED_spectrum);

if denominator_Reff == 0
    error('Sum of normalized LED intensities is zero. Check LED spectrum data.');
end
R_eff_AW = numerator_Reff / denominator_Reff; % A/W
fprintf('Effective Responsivity (R_eff): %.4f A/W\n', R_eff_AW);

% Calculate Total Incident Power on Photodiode (P_pd)
if R_eff_AW == 0
    P_pd_W = 0;
    warning('Effective responsivity is zero. This might be due to no overlap between LED spectrum and photodiode responsivity, or zero responsivity in the overlap range.');
else
    P_pd_W = Ip_A / R_eff_AW; % Watts
end
fprintf('Total Incident Power on Photodiode (P_pd): %.4e W (%.4f uW)\n', P_pd_W, P_pd_W * 1e6);
P_true = P_pd_W;

% Calculate Power per Bin (P_bin)
% Denominator is the same as for R_eff (sum of normalized LED intensities)
if denominator_Reff == 0 % Should have been caught earlier
    P_bin_W = zeros(size(LED_wavelengths));
else
    P_bin_W = P_pd_W * (LED_spectrum / denominator_Reff); % Watts per bin
end

power_per_bin = P_bin_W;

% --- 3. Display Results ---
fprintf('\n--- Power per Bin (Incident on Photodiode) ---\n');
fprintf('Wavelength (nm) | Power in Bin (uW)\n');
% fprintf('---------------------------------------\n');
% for i = 1:length(LED_wavelengths)
%     fprintf('%15.1f | %18.4f\n', LED_wavelengths(i), P_bin_W(i) * 1e6);
% end

% Sanity check: Sum of P_bin should equal P_pd
sum_P_bin_W = sum(P_bin_W);
fprintf('---------------------------------------\n');
fprintf('Sum of Power in Bins: %.4e W (%.4f uW)\n', sum_P_bin_W, sum_P_bin_W * 1e6);
fprintf('Calculated Total P_pd:  %.4e W (%.4f uW)\n', P_pd_W, P_pd_W * 1e6);
if abs(sum_P_bin_W - P_pd_W) / P_pd_W < 1e-6 % Check for small relative difference
    fprintf('Sanity check passed: Sum of P_bin is close to P_pd.\n');
else
    warning('Sanity check failed: Sum of P_bin significantly differs from P_pd. Check calculations or data.');
end

wavelengths = LED_wavelengths;