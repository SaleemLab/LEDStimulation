# Gamma Correction LUT Pipeline

Because LED forward current, MOSFET switching dynamics, and optical diffusion exhibit non-linearities (especially near switching thresholds), the system uses Look-Up Tables (LUTs) stored in microcontroller Flash memory (`PROGMEM`) to achieve true linear luminance scaling.

---

## Processing Pipeline

```mermaid
graph TD
    RawData[Raw Duty Cycle vs Measured Power Data] --> Script[Calibration/ganzfeldCalibration_v5.m]
    Script --> FitCurve[Non-Linear Curve Fitting & Monotonic Inversion]
    FitCurve --> MapScript[Calibration/mapLUTs.m]
    MapScript --> GenArray[generateGammaCorrectionLUT.m]
    GenArray --> CHeader[PROGMEM const uint16_t gammaLUT[] Array]
```

---

## 1. MATLAB Scripts & Roles

Located in [`Calibration/`](file:///d:/Code/LEDStimulation/Calibration/):

| Script | Purpose |
| :--- | :--- |
| [`ganzfeldCalibration_v5.m`](file:///d:/Code/LEDStimulation/Calibration/ganzfeldCalibration_v5.m) | Main calibration script: interfaces with serial port, coordinates `gc` sweeps, reads power meter data, and computes luminance curves. |
| [`mapLUTs.m`](file:///d:/Code/LEDStimulation/Calibration/mapLUTs.m) | Maps measured power outputs to linearized 8-bit or 16-bit register values; performs monotonic interpolation. |
| [`generateGammaCorrectionLUT.m`](file:///d:/Code/LEDStimulation/Calibration/generateGammaCorrectionLUT.m) | Converts calibrated inverse functions into formatted C header array blocks. |
| [`generateDefaultGammaCorrectionLUT.m`](file:///d:/Code/LEDStimulation/Calibration/generateDefaultGammaCorrectionLUT.m) | Generates idealized linear/identity fallback tables. |

---

## 2. Generating & Deploying a New LUT

1. **Collect Sweep Data:** Run `ganzfeldCalibration_v5.m` to measure the luminance response curves for Channel A (e.g. Green) and Channel B (e.g. UV).
2. **Compute Inverse Mapping:**
   ```matlab
   % In MATLAB
   lut = mapLUTs(measured_power, duty_cycles);
   ```
3. **Export C Code Array:**
   ```matlab
   generateGammaCorrectionLUT(lut, 'ChA_GammaLUT');
   ```
4. **Embed in Firmware:**
   - Copy the generated C array into the firmware header/sketch (e.g., `StimulusControlViaSerial_Leonardo_v8.ino`).
   - Flash the updated firmware to the microcontroller.
5. **Runtime Selection:** Use `agc, 1` or `agc, 2` to switch between calibrated LUT profiles for different ambient/attenuation modes.
