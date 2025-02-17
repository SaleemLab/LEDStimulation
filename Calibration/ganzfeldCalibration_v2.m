%% Ganzfeld calibration v2

%% define constants
lambda_s = 360; % UV-sensitive S-opsin cones
lambda_m = 510; % GREEN-senstiive M-opsin cones
h        = 4.135667E-15;    % Planck's constant [eV*s]
c        = 299792458;       % speed of light [m/s]
eV_per_J = 6.242E+18;       % [eV] per [J]
ac_um2   = 0.2;             % cone OS light collection area [µm^2], see [Nikonov et al., 2006] (http://www.ncbi.nlm.nih.gov/pubmed/16567464) for details. This is an experimentally determined value, e.g. for wt mouse cones that is fully dark-adapted, a value of 0.2 can be assumed.
ar_um2   = 0.5;             % rod OS light collection area, see above. A value of 0.5 is considered realistic.

%% Define photoreceptors

mouse_opsins = readtable('mouse_cone_opsins.txt','ReadRowNames',false);
mouse_opsins.Properties.VariableNames = {'Sopsin', 'Mopsin', 'lambda'};


%% load mouse lens transmission

mouse_transmission = readtable('Transmission_mouse_eye.txt');
mouse_transmission.Properties.VariableNames = {'rel_transmission','lambda'};
mouse_transmission.rel_transmission = mouse_transmission.rel_transmission/100;

figure
plot(mouse_transmission.lambda, mouse_transmission.rel_transmission,'k')

xlabel('wavelength')
ylabel('Relative transmission')


%% load LED spectrum from datasheets (grabit) and do some basic processing

load('UV_spec.mat')
load('green_spec.mat')
delta_lambda=1;

UVspec = UV_LED_spectrum_grabit; clear UV_LED_spectrum_grabit
greenspec = green_LED_spectrum_grabit; clear green_LED_spectrum_grabit

% sort grabit points in asending order of wavelength, then resample
[UV_wavelength, UV_sortidx] = sort(UVspec(:,1));
UV_power = UVspec(UV_sortidx,2);
UV_wavelength_resample = [floor(min(UV_wavelength)):delta_lambda:ceil(max(UV_wavelength))];
UV_power_resample = interp1(UV_wavelength, UV_power, UV_wavelength_resample, 'spline'); % Y

% do the same as above for green
[green_wavelength, green_sortidx] = sort(greenspec(:,1));
green_power = greenspec(green_sortidx,2);
green_wavelength_resample = [floor(min(green_wavelength)):delta_lambda:ceil(max(green_wavelength))];
green_power_resample = interp1(green_wavelength, green_power, green_wavelength_resample, 'spline'); % Y

%% process power meter readings to scale LED spectra - UV

lambda_meas = 385; % specified measurement wavelength for power meter
P_total = 1e-6; % output of power meter in W

% get LED spectrum and normalise (not strictly necesary to normalise)
LED_wavelengths=UV_wavelength_resample;
LED_spectrum=UV_power_resample;
% Normalize the LED spectrum
% Calculate the area under the LED spectrum using trapz
area_LED = trapz(LED_wavelengths, LED_spectrum);
LED_spectrum_norm = LED_spectrum / area_LED;  % Now, trapz(wavelengths, normSpec) = 1


% get calibration data for the photodiode
P_meter_cal = readtable('S120VC_Responsivity.xlsx');
P_meter_cal_wavelengths = P_meter_cal.Wavelength_nm_;
P_meter_cal_responsivity = P_meter_cal.Responsivity_mA_W_;

[UV_wavelengths, UV_power, UV_P_true, UV_correction_factor] = getLEDSpectraFromPowerMeter(...
    P_total, lambda_meas, LED_wavelengths, LED_spectrum, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity);

%% process power meter readings to scale LED spectra - GREEN

lambda_meas = 525; % specified measurement wavelength for power meter
P_total = 1e-6; % output of power meter in W

% get LED spectrum and normalise (not strictly necesary to normalise)
LED_wavelengths=green_wavelength_resample;
LED_spectrum=green_power_resample;
% Normalize the LED spectrum
% Calculate the area under the LED spectrum using trapz
area_LED = trapz(LED_wavelengths, LED_spectrum);
LED_spectrum_norm = LED_spectrum / area_LED;  % Now, trapz(wavelengths, normSpec) = 1


% get calibration data for the photodiode
P_meter_cal = readtable('S120VC_Responsivity.xlsx');
P_meter_cal_wavelengths = P_meter_cal.Wavelength_nm_;
P_meter_cal_responsivity = P_meter_cal.Responsivity_mA_W_;

[green_wavelengths, green_power, green_P_true, green_correction_factor] = getLEDSpectraFromPowerMeter(...
    P_total, lambda_meas, LED_wavelengths, LED_spectrum, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity);

%% zero-pad the LED emissions to

full_range = 300:699;
% pad green
idx = find(~ismember(full_range, green_wavelengths));
zero_pad = zeros(size(idx));

green_wavelengths_pad = [green_wavelengths, full_range(idx)];
green_power_pad = [green_power, zero_pad];
[~, sort_idx] = sort(green_wavelengths_pad);
green_wavelengths_pad = green_wavelengths_pad(sort_idx);
green_power_pad = green_power_pad(sort_idx);

% pad uv
idx = find(~ismember(full_range, UV_wavelengths));
zero_pad = zeros(size(idx));

UV_wavelengths_pad = [UV_wavelengths, full_range(idx)];
UV_power_pad = [UV_power, zero_pad];
[~, sort_idx] = sort(UV_wavelengths_pad);
UV_wavelengths_pad = UV_wavelengths_pad(sort_idx);
UV_power_pad = UV_power_pad(sort_idx);

%% get power at pupil
% correct for difference in size of pupil vs photodiode colleciton area

pupil_area = 0.21;  %stationary, mm^2
pupil_area = 1.91;  %running

sensor_area_mm2 = 9.7;
sensor_pupil_size_correction = pupil_area/sensor_area_mm2;

green_power_eye = green_power_pad*sensor_pupil_size_correction;
UV_power_eye = UV_power_pad*sensor_pupil_size_correction;

%% get light power at the retina 

% first account for filtering by the lens

green_power_lens = green_power_eye.*mouse_transmission.rel_transmission;
UV_power_lens = UV_power_eye.*mouse_transmission.rel_transmission;

% account for ratio of pupil area to retinal area
eye_axial_len_mm = 3;
ret_area = 0.6 *(eye_axial_len_mm/2)^2 * pi *4;
pupil_to_retina  = pupil_area/ret_area;

green_power_retina = green_power_eye*pupil_to_retina;
UV_power_retina = UV_power_eye*pupil_to_retina;

figure, hold on
plot(full_range, green_power_retina,'g')
yyaxis left
plot(full_range, UV_power_retina,'c')
xlabel('Wavelength (nm)')
ylabel('Power (W)')
title('Light hitting the retina')

yyaxis right
% mouse opsins
plot(mouse_opsins.lambda, mouse_opsins.Sopsin, 'b--')
plot(mouse_opsins.lambda, mouse_opsins.Mopsin, 'g--')

%% map electrical power to photon flux

h        = 4.135667E-15;    % Planck's constant [eV*s]
c        = 299792458;       % speed of light [m/s]
eV_per_J = 6.242E+18;       % [eV] per [J]
ac_um2   = 0.2;             % cone OS light collection area [µm^2], see [Nikonov et al., 2006] (http://www.ncbi.nlm.nih.gov/pubmed/16567464) for details. This is an experimentally determined value, e.g. for wt mouse cones that is fully dark-adapted, a value of 0.2 can be assumed.
ac_m2 = ac_um2*1E-12;
ar_um2   = 0.5;             % rod OS light collection area, see above. A value of 0.5 is considered realistic.
wavelengths = full_range;
wavelengths_m = full_range*1E-9;


P_phi = (green_power_retina.*wavelengths_m*eV_per_J)/(c*h);

photoIso_per_nm = P_phi(:) .* ac_m2 .* mouse_opsins.Mopsin(:)

R_star = trapz(wavelengths_nm, photoIso_per_nm)
% %% convert power of LEDs to photon flux density
% 
% h = 6.626e-34;        % Planck's constant in J*s
% c        = 299792458;       % speed of light [m/s]
% wavelengths = full_range(:);
% wavelengths_m = wavelengths_nm * 1e-9;
% 
% % Calculate photon flux (photons/s/nm) at each wavelength
% photonFlux = (green_power_retina(:) .* wavelengths_m(:)) / (h * c);
% 
% plot(full_range, photonFlux)
% 
% %% account for photoreceptors collection area
% 
% A_eff = 2e-13;
% photonFlux_captured = photonFlux * A_eff;
% 
% photoisomerisations_per_nm = photonFlux_captured(:).*mouse_opsins.Mopsin(:);
% 
% R_star_total = trapz(wavelengths_nm, photoisomerizations_per_nm')