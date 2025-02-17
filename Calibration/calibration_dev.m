%% LED calibration for mouse in vivo recordings
% based on https://github.com/eulerlab/open-visual-stimulator/blob/master/calibration_mouse_invivo/stimulator_calibration_in_vivo.ipynb

% 1. We first map electrical power (P_el in [W]) to photon flux (P_phi in
% [photons/s]):

% p_phi(lambda) = { (P_el(lambda) * a * lambda * 10^-9) / (c*h) } * { 1 / mu_lens2cam(lambda) }

% where lambda is peak wavelength of photoreceptor spectral sensitivity

% 2. Then convert photon flux to photoisomerisation rate (R_iso in [P*/cone/s])

% R_iso(lambda) = { P_phi(lambda) / A_stim } * A_collect * T(lambda) * Rpup2ret

% where A_stim is the area that is illuminated on the power meter and
% A_collect = 0.2um^2, the photoreceptors OS light collection area

%% define constants
lambda_s = 360; % UV-sensitive S-opsin cones
lambda_m = 510; % GREEN-senstiive M-opsin cones
h        = 4.135667E-15;    % Planck's constant [eV*s]
c        = 299792458;       % speed of light [m/s]
eV_per_J = 6.242E+18;       % [eV] per [J]
ac_um2   = 0.2;             % cone OS light collection area [µm^2], see [Nikonov et al., 2006] (http://www.ncbi.nlm.nih.gov/pubmed/16567464) for details. This is an experimentally determined value, e.g. for wt mouse cones that is fully dark-adapted, a value of 0.2 can be assumed.
ar_um2   = 0.5;             % rod OS light collection area, see above. A value of 0.5 is considered realistic.
A_detect_um2 = 9700*9700;            % Thorlabs photodiode active detector area

%% Define photoreceptors

% load data

mouse_opsins = readtable('mouse_cone_opsins.txt','ReadRowNames',false);
mouse_opsins.Properties.VariableNames = {'Sopsin', 'Mopsin', 'lambda'};


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


%% plot normalised photoreceptor sensitivity and led spectra
figure, hold on
plot(UV_wavelength_resample, UV_power_resample/100, 'b')
plot(green_wavelength_resample, green_power_resample/100, 'g')
plot(mouse_opsins.lambda, mouse_opsins.Sopsin, 'b--')
plot(mouse_opsins.lambda, mouse_opsins.Mopsin, 'g--')

legend({'UV LED', 'GREEN LED', 'Sopsin', 'Mopsin'})


%% process power meter readings to scale LED spectra - UV

lambda_meas = 385; % specified measurement wavelength for power meter
P_total = 100; % output of power meter in nW

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
P_total = 100; % output of power meter in nW

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

%% correct for pupil size vs power meter size
% dont get why this would be necessary? correct for sensor : pupil size
% later
% Sensor active area: 9.7 mm^2 converted to cm^2
%sensor_area_mm2 = 9.7;
%sensor_area_cm2 = sensor_area_mm2 / 100;  % 0.097 cm^2

% For a square sensor, compute an equivalent circular radius (in cm)
% r_eff = sqrt(sensor_area_cm2 / pi);  % ~0.176 cm
% 
% % Estimate L from the technical drawing (distance from front aperture to active area)
% % Based on the S120VC drawing, L is roughly 1 mm = 0.1 cm.
% L = 0.1;  % in cm
% 
% % Calculate the half-angle (in radians) subtended by the sensor
% theta = atan(r_eff / L);  
% theta_deg = theta * (180/pi);
% 
% % Distance from sensor (or pinhole) to the screen (in cm)
% d = 0;  % example value, adjust as needed
% 
% % Compute the projected radius on the screen:
% r_screen = tan(theta) * d;
% 
% % Compute the effective (projected) area on the screen (assume circular footprint)
% A_proj = pi * r_screen^2;
% 
% % Screen Dimensions
% screen_x = 20;  % cm
% screen_y = 20;  % cm
% A_screen = screen_x * screen_y;


% coverage_correction = A_screen / sensor_area_cm2

%% plot photoreceptor sensitivity with LED power emissison

figure, hold on

% mouse opsins
yyaxis left
plot(mouse_opsins.lambda, mouse_opsins.Sopsin, 'b--')
plot(mouse_opsins.lambda, mouse_opsins.Mopsin, 'g--')

% corrected LED emission spectrum
yyaxis right
plot(UV_wavelengths, UV_power, 'b-')
plot(green_wavelengths, green_power, 'g-')

%% load mouse lens transmission

mouse_transmission = readtable('Transmission_mouse_eye.txt');
mouse_transmission.Properties.VariableNames = {'rel_transmission','lambda'};

figure
plot(mouse_transmission.lambda, mouse_transmission.rel_transmission,'k')
xlabel('wavelength')
ylabel('Relative transmission')
uv_transmission = mouse_transmission.rel_transmission(find(mouse_transmission.lambda==360));
green_transmission = mouse_transmission.rel_transmission(find(mouse_transmission.lambda==510));



% lens_transmission_ratio = green_transmission/uv_transmission;

%% calculate Rpup2ret (stationary or locomoting)

eye_axial_len_mm = 3;
ret_area = 0.6 *(eye_axial_len_mm/2)^2 * pi *4;

pupil_area = 0.21;  %stationary
pupil_area = 1.91;  %running

pupil_to_retina = pupil_area/ret_area;

%% correct for pupil sizer compare to sensor size

sensor_area_mm2 = 9.7;
sensor_pupil_size_correction = pupil_area/sensor_area_mm2;


%% correct LED emission spectrum based on eye transmission and pupil size

% green
idx = find(ismember(mouse_transmission.lambda,green_wavelengths(:)));
green_spec_corr = green_power'.*mouse_transmission.rel_transmission(idx);
green_spec_corr = green_spec_corr/100;
green_spec_corr = green_spec_corr*pupil_to_retina*sensor_pupil_size_correction;

% uv
idx = find(ismember(mouse_transmission.lambda,UV_wavelengths(:)));
uv_spec_corr = UV_power'.*mouse_transmission.rel_transmission(idx);
uv_spec_corr = uv_spec_corr/100;
uv_spec_corr = uv_spec_corr*pupil_to_retina*sensor_pupil_size_correction;

%% add 0s to pad spectrum to full range (300-700nm)

full_range = 300:699;
% pad green
idx = find(~ismember(full_range, green_wavelengths));
zero_pad = zeros(size(idx));

green_wavelengths_pad = [green_wavelengths, full_range(idx)];
green_power_pad = [green_spec_corr', zero_pad];
[~, sort_idx] = sort(green_wavelengths_pad);
green_wavelengths_pad = green_wavelengths_pad(sort_idx);
green_power_pad = green_power_pad(sort_idx);

% pad uv
idx = find(~ismember(full_range, UV_wavelengths));
zero_pad = zeros(size(idx));

UV_wavelengths_pad = [UV_wavelengths, full_range(idx)];
UV_power_pad = [uv_spec_corr', zero_pad];
[~, sort_idx] = sort(UV_wavelengths_pad);
UV_wavelengths_pad = UV_wavelengths_pad(sort_idx);
UV_power_pad = UV_power_pad(sort_idx);

%% plot photoreceptor sensitivity with corrected LED emission power

figure, hold on

% mouse opsins
yyaxis left
plot(mouse_opsins.lambda, mouse_opsins.Sopsin, 'b--')
plot(mouse_opsins.lambda, mouse_opsins.Mopsin, 'g--')
ylabel('Photoreceptor rel. sensitivity')
% corrected LED emission spectrum
yyaxis right
plot(UV_wavelengths_pad, UV_power_pad, 'b-')
plot(green_wavelengths_pad, green_power_pad, 'g-')
ylabel('Eye-corrected LED power (nW)')

%% Determine effective photoreceptor stimulation

% get relative excitation for each photoreceptor/LED combination
for iPR = 1:2
    for iLED = 1:2
        if iPR==1, PR = mouse_opsins.Sopsin(:); end
        if iPR==2, PR = mouse_opsins.Mopsin(:); end
        if iLED==1, LED_spec_corr = UV_power_pad(:); end
        if iLED==2, LED_spec_corr = green_power_pad(:); end

        tempSpect = PR.*LED_spec_corr;
        A_overlap = trapz(tempSpect);
        A_LED = trapz(LED_spec_corr);
        rel_exc(iPR, iLED) = A_overlap/A_LED;
    end
end

UV_green_power_ratio = rel_exc(2,2)/rel_exc(1,1)

%% get photoisomerisation rates

wavelength = 300:699; wavelength=wavelength(:);
A_detect_um2 = 9.7e+6;

for iLED = 1:2

    if iLED==1, LED_spec_corr = UV_power_pad(:); end
    if iLED==2, LED_spec_corr = green_power_pad(:); end

    % Convert energy flux from [nW] (=readout of spectrometer) into [eV/s]        pow_eflux = LED_spec_corr *1E-9 *eV_per_J;
    pow_eflux = LED_spec_corr * 1E-9 * eV_per_J;

    % Calculate the wavelength-dependent photon energy `Q` in [eV]
    pow_Q = (c*h)./(wavelength*1E-9);

    % Divide energy flux by the photon energy to get the photon flux `phi`[photons/s]
    % and then photon flux density `E` [photons/s /µm^2]
    pow_phi = pow_eflux./pow_Q;
    pow_E = pow_phi ./ A_detect_um2;


    for iPR = 1:2
        if iPR==1, PR = mouse_opsins.Sopsin(:); end
        if iPR==2, PR = mouse_opsins.Mopsin(:); end

        A_collect = ac_um2;
        photon_rate = pow_E*A_collect;
        photoiso_rate = photon_rate*rel_exc(iPR,iLED);
        photoiso_rate_total = sum(photoiso_rate)*1E-9;

        photoiso(iLED,iPR).photoiso_rate = photoiso_rate;
        photoiso(iLED,iPR).photoiso_rate_total = photoiso_rate_total;



    end



end




%% calculate LED / opsin transmission

% uv_idx = find(mouse_opsins.lambda==385); %peak Sopsin lambda
% sOpsin_UVLED_transmission = mouse_opsins.Sopsin(uv_idx);
%
% green_idx = find(mouse_opsins.lambda==510); %peak Sopsin lambda
% MOpsin_greenLED_transmission = mouse_opsins.Mopsin(green_idx);
%
% led_transmission_ratio = MOpsin_greenLED_transmission/sOpsin_UVLED_transmission;
%
% %% calulate required power ratio
%
% UV_power_multiplier = lens_transmission_ratio * led_transmission_ratio;






%% From Nauhaus lab 2021 paper,https://www.nature.com/articles/s41598-021-90650-4#Sec16

% to estimate photoiso rates of rods and cones, convert radiance from the

% I(lambda) = T(lambda) * Rad(lambda) * pupil:retina ratio
% where T(lambda) is transmission, Rad(lambda) is radiance of the display 
% which is power / Area * steradians

% convert to quanta
% Iq(lambda) = I(lambda)*lambda/(c*h)
% photoIsoimerisation, take dot product between Iq(lambda) and absorption
% spectrum and scale by delta(lambda) - sample period of specturm (we use
% 1)
%R*/receptor/s = delta(lambda)


% A = 0.0097^2; % collection area of photodetector
% d= 0.17; % distance between eye and display, in cm
% sr = A/d^2;
% power = 0.001; % in Watts
% 
% % convert UV power to radiance
% L = UV_power_pad*10^-9/(A*sr);
% 
% 
% % retinal irradiance
% I = L'.*(mouse_transmission.rel_transmission/100)*pupil_to_retina;
% 
% Iq = I' /powQ


%%
% Constants
h = 6.626e-34;      % Planck's constant [J*s]
c = 3.0e8;          % Speed of light [m/s]

% Example inputs:
% wavelengths in nm (vector)
wavelengths_nm = green_wavelengths_pad(:); %linspace(400, 700, 301);  % 400 nm to 700 nm in 1 nm steps
% Convert to meters
wavelengths_m = wavelengths_nm * 1e-9;

% Spectral power P (in Watts per nm) reaching the retina at each wavelength
% (Replace this with your actual measured/simulated data)
P_lambda = green_power_pad(:);  % should be a vector of same length as wavelengths_nm

% Photoreceptor parameters
A_eff = 2e-13; %1e-12;   % Effective collecting area in m^2 (example value; adjust as needed)
% Spectral sensitivity S(lambda): unitless weighting factor (vector, same length)
S_lambda = mouse_opsins.Mopsin(:);  % e.g., from opsin absorption data

% Calculate photon flux at each wavelength (photons/s per nm)
photonFlux = (P_lambda .* wavelengths_m) ./ (h * c);

% Calculate photoisomerizations per nm
photoisomerizations_per_nm = photonFlux .* A_eff .* S_lambda;

% Integrate over the wavelength range to get total photoisomerizations (per second)
% Use trapz for numerical integration; wavelengths in meters
R_star = trapz(wavelengths_m, photoisomerizations_per_nm);
R_star2 = sum(photoisomerizations_per_nm)
fprintf('Estimated photoisomerizations: %.2f per second\n', R_star);
