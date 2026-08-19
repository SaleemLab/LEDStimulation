# Optical Calibration Guide

Precise visual neuroscience and psychophysics require quantitative radiometric and photometric calibration. This guide covers the procedure for measuring absolute power, determining baseline operating points, and linearizing stimulus output.

---

## Calibration Workflow

```mermaid
graph LR
    Setup[Mount Sensor in Enclosure] --> PotAlign[Potentiometer Baseline Adjustment]
    PotAlign --> GammaSweep[Run Automated 'gc' Sweep]
    GammaSweep --> MATLAB[Process Data in MATLAB]
    MATLAB --> GenLUT[Generate C LUT Array]
    GenLUT --> UploadFW[Embed LUT in Firmware & Flash]
```

---

## 1. Required Equipment

- **Optical Power Meter Console:** (e.g. Thorlabs PM100D or PM400).
- **Calibrated Photodiode Sensor:**
  - Must have appropriate spectral sensitivity covering the UV (365–405 nm) and Visible (500–550 nm) ranges.
  - Sensor calibration file / wavelength responsiveness table loaded on the console.
- **Mounting:** Position the sensor at the exact plane of the subject's eye / retina relative to the diffusion surface.

---

## 2. Potentiometer & Baseline Setup

1. **Hardware Mode:** The driver allows peak current limiting and analog gate scaling via onboard potentiometers.
2. **Streaming Sensor Feedback:**
   - Send the `ana` command via Serial Monitor or Bonsai to stream continuous readings from pins `A0` and `A1`.
   - Adjust the multi-turn potentiometers until the desired baseline voltage / peak current is reached.
   - Send `done` to terminate the stream.
3. Record the potentiometer positions in [`PotentiometerSettingsNotes.txt`](file:///d:/Code/LEDStimulation/PotentiometerSettingsNotes.txt).

---

## 3. Automated Gamma Measurement Sweep

1. Open a serial terminal or run the MATLAB calibration runner ([`Calibration/ganzfeldCalibration_v5.m`](file:///d:/Code/LEDStimulation/Calibration/ganzfeldCalibration_v5.m)).
2. Trigger the automated step sequence using the `gc` command:
   ```text
   gc, 2, 250, 3
   ```
   *(Steps duty cycle in 2% increments from 0% to 100%, pausing 250ms per step, averaging over 3 repetitions).*
3. Record the synchronized optical power measurements ($\mu\text{W}$) for each channel independently.
