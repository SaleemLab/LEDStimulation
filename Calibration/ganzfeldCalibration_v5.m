%% ganzfeld calibration v3


%% select measurement device

measurementDevice = 'PDA100A2'; % 'PDA100A2, S130VC, S120VC
load('PDA100A2_gainTable.mat')

%% Define constants

% h: Planck's constant [eV*s]
% c: speed of light [m/s]
% eV_per_J: conversion factor ([eV] per [J])
% ac_um2: cone OS light collection area [µm^2], see [Nikonov et al., 2006] (http://www.ncbi.nlm.nih.gov/pubmed/16567464) for details. This is an experimentally determined value, e.g. for wt mouse cones that is fully dark-adapted, a value of 0.2 can be assumed.
% ar_um2: rod OS light collection area, see above. A value of 0.5 is considered realistic.

h        = 4.135667E-15; % Planck's constant [eV*s]
c        = 299792458;    % speed of light [m/s]
eV_per_J = 6.242E+18;    % [eV] per [J]
ac_um2   = 0.2;     
ar_um2   = 0.5;

%% Define photoreceptors

mouse_opsins = readtable('mouse_cone_opsins.txt','ReadRowNames',false);
mouse_opsins.Properties.VariableNames = {'Sopsin', 'Mopsin', 'lambda'};

% M-cone
PRs(1).name = "mouse_M_cone";
PRs(1).peak_nm = 511;
PRs(1).collecArea_um2 = ac_um2;
PRs(1).spect = mouse_opsins.Mopsin;
PRs(1).col = 'g';
% S-cone
PRs(2).name = "mouse_S_cone";
PRs(2).peak_nm = 360;
PRs(2).collecArea_um2 = ac_um2;
PRs(2).spect = mouse_opsins.Sopsin;
PRs(2).col = 'm';
% Rod
PRs(3).name = "mouse_rod";
PRs(3).peak_nm = 510;
PRs(3).collecArea_um2 = ar_um2;
PRs(3).spect = PRs(1).spect;
PRs(3).col = 'r';


%% Define light sources

%load('UV_spec.mat')
%load('green_spec.mat')
load('D:\Code\LEDStimulation\Calibration\LEDdata\GREEN_grabit.mat')
load('D:\Code\LEDStimulation\Calibration\LEDdata\UV_grabit.mat')

delta_lambda=1;

UVspec = UV_grabit; clear UV_grabit
greenspec = GREEN_grabit; clear green_grabit

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

% zero-pad both spectra to full wavelength range
wavelengths_nm = 300:699;

% pad uv
idx = find(~ismember(wavelengths_nm, UV_wavelength_resample));
zero_pad = zeros(size(idx));

UV_wavelengths_pad = [UV_wavelength_resample, wavelengths_nm(idx)];
UV_power_pad = [UV_power_resample, zero_pad];
[~, sort_idx] = sort(UV_wavelengths_pad);
UV_wavelengths_pad = UV_wavelengths_pad(sort_idx);
UV_power_pad = UV_power_pad(sort_idx);

% pad green
idx = find(~ismember(wavelengths_nm, green_wavelength_resample));
zero_pad = zeros(size(idx));

green_wavelengths_pad = [green_wavelength_resample, wavelengths_nm(idx)];
green_power_pad = [green_power_resample, zero_pad];
[~, sort_idx] = sort(green_wavelengths_pad);
green_wavelengths_pad = green_wavelengths_pad(sort_idx);
green_power_pad = green_power_pad(sort_idx);

% name: name of LED/filter combinations, used later for plots etc.
% peak_nm: peak wavelength of LED/filter combination in [nm]
% LED_spect: spectrum of LED (same x range as the opsin spectra, that is from 300 to 699 nm, and with 1-nm resolution).
% spect_nw: measured LED/filter spectrum in [nW]
% spect_nw_norm: spect_nw peak-normalized

LEDs(1).name = "green";
LEDs(1).peak_nm = wavelengths_nm(green_power_pad==max(green_power_pad));
LEDs(1).target_PR = "mouse_M_cone";
LEDs(1).spect_nw = nan;
LEDs(1).spect_nw_norm = green_power_pad./max(green_power_pad);

LEDs(2).name = "UV";
LEDs(2).peak_nm = wavelengths_nm(UV_power_pad==max(UV_power_pad));
LEDs(2).target_PR = "mouse_S_cone";
LEDs(2).spect_nw = nan;
LEDs(2).spect_nw_norm = UV_power_pad./max(UV_power_pad);

%% plot photoreceptor sensitivity and LED spectra

figure, hold on
for iPR=1:2
    plot(wavelengths_nm, PRs(iPR).spect,'color', PRs(iPR).col);
    plot(wavelengths_nm, LEDs(iPR).spect_nw_norm, 'Color', PRs(iPR).col, 'LineStyle', ':')
end
xlabel('wavelength (nm)'), ylabel('Rel. sensitivity'), 
legend({PRs(1).name, LEDs(1).name, PRs(2).name, LEDs(2).name});
defaultAxesProperties(gca, false)


%% load appropriate responsivity data and detection area of photodiode sensor
switch measurementDevice
    case 'S130VC'
        P_meter_cal = readtable(fullfile('D:\Code\LEDStimulation\Calibration\', 'PowerMeterCalibrationFiles', 'S130VC_Responsivity.xlsx'));
        P_meter_cal_wavelengths = P_meter_cal.Wavelength_nm_;
        P_meter_cal_responsivity = P_meter_cal.Responsivity_mA_W_; % without filter

        A_detect_um2 = 9700^2; %Thor labs photodiode active receptor area

    case 'S120VC'
        P_meter_cal = readtable(fullfile('D:\Code\LEDStimulation\Calibration\', 'PowerMeterCalibrationFiles', 'S120VC_Responsivity.xlsx'));
        P_meter_cal_wavelengths = P_meter_cal.Wavelength_nm_;
        P_meter_cal_responsivity = P_meter_cal.Responsivity_mA_W_; % without filter

        A_detect_um2 = 9700^2; %Thor labs photodiode active receptor area

    case 'S120C'
        P_meter_cal = readtable(fullfile('D:\Code\LEDStimulation\Calibration\', 'PowerMeterCalibrationFiles', 'S120C_Responsivity.xlsx'));
        P_meter_cal_wavelengths = P_meter_cal.Wavelength_nm_;
        P_meter_cal_responsivity = P_meter_cal.Responsivity_mA_W_; % without filter

        A_detect_um2 = 9700^2; %Thor labs photodiode active receptor area

    case 'PDA100A2'

        P_meter_cal = readtable(fullfile('D:\Code\LEDStimulation\Calibration\', 'PowerMeterCalibrationFiles', 'PDA100A2_Responsivity.xlsx'));
        P_meter_cal_wavelengths = P_meter_cal.Wavelength_nm_;
        P_meter_cal_responsivity = P_meter_cal.Responsivity_A_W_; % without filter

        A_detect_um2 = 75400000; % 75.4mm^2
        
end

%% get correct emission spectra of LEDs based on responsivity of measurment device


if strcmp(measurementDevice, 'PDA100A2') % thorlabs photodiode

    % green measurements
    V_measured = 7.222; %*10^-1; ^-2
    V_dark=0;
    gainVal = PDA100A2_gainTable.("Hi-Z gain (V/A)")(PDA100A2_gainTable.("dB Setting")==70);


    [green_wavelengths, green_power, green_P_true] = getLEDSpectraFromPhotodiode(...
    V_measured, V_dark, wavelengths_nm,  LEDs(1).spect_nw_norm, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity, gainVal);
    
    % uv measurements
    V_measured = 9.10; % *10^-1; ^-2*10^-1;
    V_dark=0;
    gainVal = PDA100A2_gainTable.("Hi-Z gain (V/A)")(PDA100A2_gainTable.("dB Setting")==70);


    [UV_wavelengths, UV_power, UV_P_true] = getLEDSpectraFromPhotodiode(...
    V_measured, V_dark, wavelengths_nm,  LEDs(2).spect_nw_norm, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity, gainVal);

else % thorlabs power meter
% use Thor-labs power meter sensor responsivity calibration data to
% convert measured power into a correctedf emission spectra for each light 
% source in nW.


% process power meter readings to scale LED spectra - GREEN

lambda_meas = 525; % specified measurement wavelength for power meter
P_total = 6.5E-6; % output of power meter in W

[wavelengths, green_power, green_P_true, green_correction_factor] = getLEDSpectraFromPowerMeter(...
    P_total, lambda_meas, wavelengths_nm, LEDs(1).spect_nw_norm, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity);

% green_correction_factor



%process power meter readings to scale LED spectra - UV

lambda_meas = 370; % specified measurement wavelength for power meter
P_total = 15E-6; % output of power meter in W

[UV_wavelengths, UV_power, UV_P_true, UV_correction_factor] = getLEDSpectraFromPowerMeter(...
    P_total, lambda_meas, wavelengths_nm, LEDs(2).spect_nw_norm, ...
    P_meter_cal_wavelengths, P_meter_cal_responsivity);

end

LEDs(1).spect_nw = green_power*1e9;
LEDs(1).measured_power = green_P_true; % reference for the spectrum

LEDs(2).spect_nw = UV_power*1e9; % convert to nw to match OVS notebook
LEDs(2).measured_power = UV_P_true; % reference for the spectrum

%% plot corrected-power spectra of LEDs 
figure, hold on
plot(wavelengths_nm, LEDs(1).spect_nw, 'Color', PRs(1).col)
plot(wavelengths_nm, LEDs(2).spect_nw, 'Color', PRs(2).col)
defaultAxesProperties(gca, false)
xlabel('wavelength (nm)'), ylabel('power (nW)')

%% correct emission spectra based on properties of mouse eye
% 
% load mouse lens transmission

mouse_transmission = readtable('Transmission_mouse_eye.txt');
mouse_transmission.Properties.VariableNames = {'rel_transmission','lambda'};
mouse_transmission.rel_transmission = mouse_transmission.rel_transmission/100;

% figure
% plot(mouse_transmission.lambda, mouse_transmission.rel_transmission,'k')
% 
% xlabel('wavelength')
% ylabel('Relative transmission')

% get pupil:retina ratio
% pupil_area = 0.21;  %stationary, mm^2
pupil_area = 1.91;  %running

% account for ratio of pupil area to retinal area
eye_axial_len_mm = 3;
ret_area = 0.6 *(eye_axial_len_mm/2)^2 * pi *4;
pupil_to_retina  = pupil_area/ret_area;

for iLED = 1:2
    LEDs(iLED).spect_nw_corr = LEDs(iLED).spect_nw .* mouse_transmission.rel_transmission' .* pupil_to_retina;
    LEDs(iLED).spect_nw_corr_norm = LEDs(iLED).spect_nw_corr ./max(LEDs(iLED).spect_nw_corr );
end


%% Determine effective photoreceptor stimulation

for iLED = 1:2
    for iPR = 1:3
        tempSpect = PRs(iPR).spect.*LEDs(iLED).spect_nw_corr_norm(:);
        A_overlap = trapz(tempSpect);
        A_LED = trapz(LEDs(iLED).spect_nw_corr_norm(:));
        LEDs(iLED).effect_onPR(iPR) = A_overlap/A_LED;
    end
end

fprintf('Relative co-excitation:\n')
for iLED = 1:2
    for iPR = 1:3

    fprintf('%s LED on %s: %.3f%% \n', LEDs(iLED).name, PRs(iPR).name,  LEDs(iLED).effect_onPR(iPR)*100)

    end
end
fprintf('\n')

%% Calculate and print photo-isomerization rates for all LED/filter and 

for iLED = 1:2
    
    % convert energy flux from nW to eV/s
    pow_eflux = LEDs(iLED).spect_nw_corr * 1E-9 * eV_per_J; 

    % calculate wavelength-dependent photon energy 'Q' in eV
    pow_Q     = (c*h)./(wavelengths_nm * 1E-9); % eh: is this correct with h in eV?
     
    % Divide energy flux by the photon energy to get the photon flux `phi`[photons/s] 
    %  and then photon flux density `E` [photons/s /µm^2]
    pow_phi   = pow_eflux./pow_Q;
    pow_E     = pow_phi /A_detect_um2;

    LEDs(iLED).pow_eflux = pow_eflux;
    LEDs(iLED).pow_Q = pow_Q;
    LEDs(iLED).pow_phi = pow_phi;
    LEDs(iLED).pow_E = pow_E;

    for iPR = 1:3
        % photon flux per photoreceptor 'photon_rate' in photons/s
        A_collect = PRs(iPR).collecArea_um2;
        photon_rate = LEDs(iLED).pow_E * A_collect;
        photoiso_rate = photon_rate * LEDs(iLED).effect_onPR(iPR);
        photoiso_rate_total = sum(photoiso_rate);
        
        LEDs(iLED).effect_on_PR(iPR).iso_total = photoiso_rate_total;
        LEDs(iLED).effect_on_PR(iPR).photoiso_rate = photoiso_rate;

    end
end

%% print photoisomerisation rates

fprintf('photoisomerization rates:\n')
for iLED = 1:2
    for iPR = 1:3

    fprintf('%s LED on %s: %.3f 10^3 photons/s\n', LEDs(iLED).name, PRs(iPR).name,  LEDs(iLED).effect_on_PR(iPR).iso_total*1e-3)

    end
end
fprintf('\n')


%% plot absorption spectra and photoisorates

figure, hold on
for iLED=1:2
    for iPR = 1:2
        x=wavelengths_nm;
        y=LEDs(iLED).effect_on_PR(iPR).photoiso_rate;
        % Fill area under the curve
        fill([x fliplr(x)], [y zeros(size(y))], PRs(iPR).col, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    end
end
ylabel('Isomerisation rate')
xlabel('wavelengths (nm)')
ax=gca; ax.YLim(1)=0;

yyaxis right
for iPR = 1:2
plot(wavelengths_nm, PRs(iPR).spect,'color', PRs(iPR).col, 'LineStyle','-');
end
ylabel('Photoreceptor sensitivity / LED emission spectra (normalised)')
    
for iLED = 1:2
plot(wavelengths_nm, LEDs(iLED).spect_nw_corr_norm, 'Color', PRs(iLED).col, 'LineStyle', '--')
end

legend({'green on M-cone', 'green on S-cone','UV on M-cone', 'UV on S-cone',...
    'M-cone sensitivity', 'S-cone sensitivity', 'Green LED spectra', 'UV LED spectra'})

ax=gca; ax.YLim(1)=0;