# Serial API Reference

The visual stimulator communicates with host software (e.g., Bonsai, MATLAB, Python) via USB Serial using comma-separated ASCII commands.

* **Baud Rate:** `115200`
* **Data Format:** `8` Data Bits | `No` Parity | `1` Stop Bit
* **Line Termination:** Carriage Return (`\r`) or Newline (`\n`)
* **Completion Handshake:** After executing any timed stimulus command, the firmware returns `-1\n` over Serial to notify the host that the presentation is finished.

---

## Command Quick Reference

| Command Type | Key | Syntax | Primary Firmware | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Dual-Frequency Sinewave** | `s` | `s, Duration, FreqA, FreqB, PhaseA, PhaseB, ContrastA, ContrastB` | `DDS_8bit_2freq` *(Canonical)* | Dual-channel sine flicker with independent frequencies & raised cosine window |
| **Contrast Envelope Sine** | `se` | `se, Duration, FreqA, FreqB, EnvFreq, MaxContrastA, MaxContrastB` | `DDS_8bit_2freq` *(Canonical)* | Sinusoidal flicker modulated by a sinusoidal contrast envelope |
| **Stepped Frequency Sweep** | `sfs` | `sfs, StartFreq, EndFreq, StepFreq, CyclesPerFreq, PhaseA, PhaseB, ContrastA, ContrastB` | `DDS_8bit_2freq` *(Canonical)* | Frequency sweep holding exact integer cycles per frequency step |
| **Continuous Chirp Sweep** | `fs` | `fs, Fmin, Fmax, SweepFactor, PhaseA, PhaseB, ContrastA, ContrastB` | `DDS_8bit_2freq` & `v8` | Exponential continuous frequency sweep from $F_{\text{min}} \to F_{\text{max}} \to F_{\text{min}}$ |
| **Timed DC Step** | `sdt` | `sdt, DutyA, DutyB, Duration` | `DDS_8bit_2freq` & `v8` | Fixed DC luminance step for a defined duration |
| **Static DC Output** | `sd` | `sd, DutyA, DutyB` | All Variants | Continuous static DC output ($0-100\%$) |
| **Gaussian White Noise** | `wn` | `wn, Duration, UpdateTime, MeanFrac, StdFrac` | `Leonardo_v8` *(Stochastic)* | Continuous pseudo-random Gaussian luminance noise |
| **Frozen White Noise** | `fwn` | `fwn, Duration, UpdateTime, Reps, Seed` | `Leonardo_v8` *(Stochastic)* | Deterministic, repeatable pseudo-random Gaussian sequence |
| **Switching Noise** | `cs` | `cs, UpdateTime, SwitchTime, Reps, Mean1, Cont1, Mean2, Cont2` | `Leonardo_v8` *(Stochastic)* | 2-state adaptation noise alternating contrast/mean distributions |
| **Gamma Calibration** | `gc` | `gc, Step, Wait, Reps` | All Variants | Step-wise luminance staircase for photodiode/power meter calibration |
| **Select LUT Bank** | `agc` | `agc, Index` | All Variants | Switch active gamma LUT bank (`1` = Linear, `2` = Calibrated) |
| **Channel Enable** | `useChA` / `useChB` | `useChA, 1/0` / `useChB, 1/0` | All Variants | Enable (`1`) or mute (`0`) LED channel output |
| **Analog Stream** | `ana` | `ana` | All Variants | Stream live potentiometer readings from A0/A1 (`done` to stop) |
| **Brightness Scale** | `bright` | `bright, Channel, Percent` | `Teensy41_DDS` | Dynamic master attenuation without rebuilding LUT |
| **Upload Gamma LUT** | `loadlut` | `loadlut, Channel, v0, v1, ...` | `Teensy41_DDS` | Upload 4096-point empirical LUT over serial at runtime |

---

## Detailed Command Specifications

### 1. Dual-Frequency Sinusoidal Flicker (`s`)

Generates single- or dual-channel sinusoidal flicker with independent frequency, phase, and contrast modulation on each channel. Includes a 100 ms raised cosine fade-in/fade-out window.

```text
s, Duration, FreqA, FreqB, PhaseA, PhaseB, ContrastA, ContrastB
```

* **`Duration`** *(long, ms)*: Total duration of the stimulus in milliseconds.
* **`FreqA` / `FreqB`** *(float, Hz)*: Temporal modulation frequency in Hertz for Channel A and Channel B (e.g. `2.0`, `10.5`).
* **`PhaseA` / `PhaseB`** *(float, 0.0 – 1.0)*: Relative starting phase offset ($1.0 = 2\pi\text{ rad} = 360^\circ$). For antiphase stimulation ($180^\circ$), set Channel A to `0.0` and Channel B to `0.5`.
* **`ContrastA` / `ContrastB`** *(float, 0.0 – 1.0)*: Michelson contrast / modulation depth ($0.0 = 0\%$, $1.0 = 100\%$).
* **Example:** `s, 5000, 4.0, 4.0, 0.0, 0.5, 0.8, 0.8` *(5-second 4 Hz antiphase flicker at 80% contrast)*.
* **Dual-Frequency Example:** `s, 10000, 5.0, 7.5, 0.0, 0.0, 1.0, 1.0` *(10-second dichoptic flicker: 5 Hz on Ch A, 7.5 Hz on Ch B)*.

---

### 2. Sinusoidal Flicker with Contrast Envelope (`se`)

Modulates high-frequency sinusoidal flicker with a secondary low-frequency sinusoidal contrast envelope. The envelope starts at absolute trough ($270^\circ = 0\text{ contrast}$) and smoothly increases to peak contrast and back.

```text
se, Duration, FreqA, FreqB, EnvFreq, MaxContrastA, MaxContrastB
```

* **`Duration`** *(long, ms)*: Total stimulus duration in milliseconds.
* **`FreqA` / `FreqB`** *(float, Hz)*: Carrier frequency in Hertz for Channel A and Channel B.
* **`EnvFreq`** *(float, Hz)*: Modulation frequency of the contrast envelope (e.g. `0.2` Hz for a 5-second envelope cycle).
* **`MaxContrastA` / `MaxContrastB`** *(float, 0.0 – 1.0)*: Peak modulation contrast reached at envelope maximum.
* **Example:** `se, 10000, 10.0, 10.0, 0.5, 0.9, 0.9` *(10-second 10 Hz flicker enveloped by a 0.5 Hz contrast cycle at 90% peak contrast)*.

---

### 3. Stepped Frequency Sweep (`sfs`)

Sweeps through a series of discrete frequencies from `StartFreq` to `EndFreq` in steps of `StepFreq`, holding each frequency for an exact integer number of full sine wave cycles (`CyclesPerFreq`).

```text
sfs, StartFreq, EndFreq, StepFreq, CyclesPerFreq, PhaseA, PhaseB, ContrastA, ContrastB
```

* **`StartFreq` / `EndFreq`** *(float, Hz)*: Initial and final frequencies in Hertz. If `StartFreq < EndFreq`, sweeps upward; if `StartFreq > EndFreq`, sweeps downward.
* **`StepFreq`** *(float, Hz)*: Frequency increment per step (automatically forced positive).
* **`CyclesPerFreq`** *(int)*: Number of complete sine cycles presented at each frequency step.
* **`PhaseA` / `PhaseB`** *(float, 0.0 – 1.0)*: Initial phase offsets.
* **`ContrastA` / `ContrastB`** *(float, 0.0 – 1.0)*: Modulation contrast depth.
* **Example:** `sfs, 1.0, 20.0, 1.0, 5, 0.0, 0.0, 0.8, 0.8` *(Sweeps from 1 Hz to 20 Hz in 1 Hz steps, presenting exactly 5 full cycles per step)*.

!!! tip "Stepped Sweep Duration Calculation"
    The exact total presentation duration for `sfs` can be calculated in MATLAB using [`calculateSteppedSweepDuration.m`](file:///d:/Code/LEDStimulation/ArduinoCode/calculateSteppedSweepDuration.m).

---

### 4. Continuous Exponential Frequency Chirp (`fs`)

Generates a continuous exponential frequency chirp starting at $F_{\text{min}}$, accelerating up to $F_{\text{max}}$, and smoothly returning to $F_{\text{min}}$.

```text
fs, Fmin, Fmax, SweepFactor, PhaseA, PhaseB, ContrastA, ContrastB
```

* **`Fmin` / `Fmax`** *(float, Hz)*: Minimum and maximum frequencies of the chirp sweep.
* **`SweepFactor`** *(float)*: Exponential rate factor $M$ from $f(t) = f_0 \exp(Mt)$.
* **`PhaseA` / `PhaseB`** *(float, 0.0 – 1.0)*: Initial phase offsets.
* **`ContrastA` / `ContrastB`** *(float, 0.0 – 1.0)*: Modulation contrast depths.
* **Example:** `fs, 1.0, 30.0, 0.5, 0.0, 0.0, 0.75, 0.75` *(Continuous logarithmic chirp between 1 Hz and 30 Hz)*.

---

### 5. Gaussian White Noise (`wn`) — *Leonardo v8 Firmware*

Generates continuous pseudo-random luminance frames drawn from a Gaussian distribution at a fixed frame rate.

```text
wn, Duration, UpdateTime, MeanFrac, StdFrac
```

* **`Duration`** *(long, ms)*: Total presentation duration in milliseconds.
* **`UpdateTime`** *(long, ms)*: Duration of each discrete frame in milliseconds (e.g. `10` ms for 100 Hz frame rate).
* **`MeanFrac`** *(float, 0.0 – 1.0)*: Mean baseline luminance level ($0.5 = 50\%$ duty cycle).
* **`StdFrac`** *(float, 0.0 – 1.0)*: Standard deviation of the Gaussian distribution (sets noise contrast).
* **Example:** `wn, 10000, 10, 0.5, 0.15` *(10-second Gaussian noise, 100 Hz refresh rate, 50% mean luminance, 15% standard deviation)*.

---

### 6. Frozen White Noise (`fwn`) — *Leonardo v8 Firmware*

Pre-computes a deterministic Gaussian noise sequence into a local buffer and repeats it identically across blocks to permit trial-averaged neural response analysis.

```text
fwn, Duration, UpdateTime, Reps, Seed
```

* **`Duration`** *(long, ms)*: Duration of a single frozen sequence in milliseconds.
* **`UpdateTime`** *(int, ms)*: Frame interval in milliseconds.
* **`Reps`** *(int)*: Number of times to loop the exact sequence.
* **`Seed`** *(int)*: PRNG seed integer to ensure identical random sequences across days/rigs.
* **Example:** `fwn, 2000, 10, 5, 42` *(Repeats identical 2-second Gaussian sequence 5 times using seed 42)*.

---

### 7. Contrast-Switching Adaptation Noise (`cs`) — *Leonardo v8 Firmware*

Periodically alternates between two distinct Gaussian distributions to probe dynamic sensory gain control and neural contrast adaptation kinetics.

```text
cs, UpdateTime, SwitchTime, Reps, Mean1, Cont1, Mean2, Cont2
```

* **`UpdateTime`** *(int, ms)*: Noise frame interval in milliseconds.
* **`SwitchTime`** *(int, ms)*: Duration in milliseconds before toggling to the alternative distribution.
* **`Reps`** *(int)*: Total number of distribution switch cycles.
* **`Mean1`, `Cont1`**: Mean ($0.0-1.0$) and Contrast/Std for State 1.
* **`Mean2`, `Cont2`**: Mean ($0.0-1.0$) and Contrast/Std for State 2.
* **Example:** `cs, 10, 3000, 10, 0.5, 0.05, 0.5, 0.25` *(Alternates every 3 seconds between 5% and 25% contrast for 10 cycles)*.

---

### 8. Timed Luminance Step (`sdt`)

Applies a static DC duty cycle for a specified duration before resetting to baseline.

```text
sdt, DutyA, DutyB, Duration
```

* **`DutyA` / `DutyB`** *(float, 0.0 – 100.0)*: Output duty cycle percentage for Channel A and Channel B.
* **`Duration`** *(long, ms)*: Step hold duration in milliseconds.
* **Example:** `sdt, 75.0, 0.0, 1000` *(1-second 75% luminance pulse on Channel A)*.

---

### 9. Static Output (`sd`)

Sets continuous DC duty cycle without automatic timeout.

```text
sd, DutyA, DutyB
```

* **`DutyA` / `DutyB`** *(float, 0.0 – 100.0)*: Output duty cycle percentage.
* **Example:** `sd, 0, 0` *(Turns off both LED channels)*.

---

### 10. Calibration & Utility Commands

#### Gamma Calibration Staircase (`gc`)

```text
gc, Step, Wait, Reps
```

Iterates through luminance levels from $0\%$ to $100\%$ in increments of `Step` percent, holding each level for `Wait` milliseconds, repeated `Reps` times. Used during optical power meter or photodiode calibration sweeps.

#### Active Gamma LUT Selection (`agc`)

```text
agc, Index
```

Switches active gamma correction LUT bank (`1` = Bank 1 / Linear, `2` = Bank 2 / Calibrated empirical LUT).

#### Channel Enable / Mute (`useChA` / `useChB`)

```text
useChA, 1
useChB, 0
```

Enables (`1`) or mutes (`0`) the selected output channel.

#### Analog Stream (`ana`)

```text
ana
```

Streams continuous analog sensor values from pins A0 and A1 over Serial at 100 ms intervals to monitor manual potentiometer adjustments. Send `done\n` to terminate streaming.

---

## 11. Teensy 4.1 Extensions

When running the [`LEDStimController_Teensy41_DDS.ino`](file:///d:/Code/LEDStimulation/ArduinoCode/LEDStimController_Teensy41_DDS/LEDStimController_Teensy41_DDS.ino) firmware, the following additional commands are supported:

### Master Brightness Scaling (`bright`)

```text
bright, Channel, Percent
```

Applies post-LUT master digital attenuation ($0-100\%$) to Channel `A` or `B` without needing to re-generate or re-upload gamma tables.

* **Example:** `bright, A, 50` *(Scales Channel A maximum brightness to 50%)*.

### Runtime Gamma LUT Upload (`loadlut`)

```text
loadlut, Channel, v0, v1, v2, ...
```

Uploads a 4096-entry 12-bit calibration curve directly over serial into RAM for immediate runtime linearization.
