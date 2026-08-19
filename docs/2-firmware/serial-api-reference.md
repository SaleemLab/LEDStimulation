# Serial API Reference

The stimulator communicates via USB Serial using comma-separated ASCII commands terminated by a carriage return (`\r`) or newline (`\n`).

* **Baud Rate:** `115200`
* **Data Bits:** `8` | **Parity:** `None` | **Stop Bits:** `1`

---

## Command Quick Reference

| Command Type | Key | Syntax | Description |
| :--- | :--- | :--- | :--- |
| **Sinusoidal Flicker** | `s` | `s, Duration, Freq, PhaseA, PhaseB, ContrastA, ContrastB` | Independent dual-channel sine wave flicker |
| **Gaussian White Noise** | `wn` | `wn, Duration, UpdateTime, MeanFrac, StdFrac` | Continuous Gaussian luminance noise |
| **Frozen White Noise** | `fwn` | `fwn, Duration, UpdateTime, Reps, Seed` | Repeatable deterministic Gaussian sequence |
| **Switching Noise** | `cs` | `cs, UpdateTime, SwitchTime, Reps, Mean1, Cont1, Mean2, Cont2` | Alternating 2-state adaptation noise |
| **Frequency Chirp** | `fs` | `fs, Fmin, Fmax, SweepFactor, PhaseA, PhaseB, ContrastA, ContrastB` | Exponential frequency sweep |
| **Timed Step** | `sdt` | `sdt, DutyA, DutyB, Duration` | Step luminance for a fixed duration |
| **Static Output** | `sd` | `sd, DutyA, DutyB` | Set constant DC duty cycle ($0 - 100\%$) |
| **Gamma Calibration** | `gc` | `gc, Step, Wait, Reps` | Step-wise luminance sweep for calibration |
| **Select LUT** | `agc` | `agc, Index` | Switch active gamma LUT bank (`1` or `2`) |
| **Channel Enable** | `useChA` / `useChB` | `useChA, 1/0` / `useChB, 1/0` | Enable (`1`) or mute (`0`) channel |
| **Analog Stream** | `ana` | `ana` | Stream live potentiometer readings (`done` to stop) |

---

## Detailed Command Specifications

### 1. Sinusoidal Flicker (`s`)
Generates single- or dual-channel sinusoidal flicker with independent phase and contrast modulation.
```text
s, Duration, Freq, PhaseA, PhaseB, ContrastA, ContrastB
```
* **`Duration`** *(int, ms)*: Total duration of the stimulus in milliseconds.
* **`Freq`** *(float, Hz)*: Temporal frequency in Hertz (e.g. `2.0`, `10.5`).
* **`PhaseA` / `PhaseB`** *(float, 0.0 – 1.0)*: Relative phase offset for Channel A and Channel B ($1.0 = 2\pi$ radians). For antiphase stimulation, set one channel to `0.0` and the other to `0.5`.
* **`ContrastA` / `ContrastB`** *(float, 0.0 – 1.0)*: Michelson contrast / modulation depth around the mean.
* **Example:** `s, 5000, 2.0, 0.0, 0.5, 0.8, 0.8` *(5-second 2 Hz antiphase flicker at 80% contrast)*.

---

### 2. Gaussian White Noise (`wn`)
Generates continuous pseudo-random luminance frames drawn from a Gaussian distribution.
```text
wn, Duration, UpdateTime, MeanFrac, StdFrac
```
* **`Duration`** *(int, ms)*: Total stimulus duration in milliseconds.
* **`UpdateTime`** *(int, ms)*: Duration of each discrete noise frame (e.g. `10` ms for 100 Hz refresh).
* **`MeanFrac`** *(float, 0.0 – 1.0)*: Mean baseline luminance fraction.
* **`StdFrac`** *(float, 0.0 – 1.0)*: Standard deviation of the Gaussian distribution (sets noise contrast).
* **Example:** `wn, 10000, 10, 0.5, 0.15` *(10-second noise, 100 Hz frame rate, 50% mean, 15% std)*.

---

### 3. Frozen White Noise (`fwn`)
Plays a deterministic, repeatable pseudo-random Gaussian sequence for trial-averaging neural responses.
```text
fwn, Duration, UpdateTime, Reps, Seed
```
* **`Duration`** *(int, ms)*: Duration of a single frozen sequence.
* **`UpdateTime`** *(int, ms)*: Frame duration in milliseconds.
* **`Reps`** *(int)*: Number of times to loop the exact sequence.
* **`Seed`** *(unsigned long)*: Integer PRNG seed.
* **Example:** `fwn, 2000, 10, 5, 42` *(Repeat identical 2-second noise 5 times using seed 42)*.

---

### 4. Switching White Noise (`cs`)
Alternates between two distinct Gaussian distributions (e.g., high vs. low contrast, or high vs. low mean) to probe dynamic sensory adaptation.
```text
cs, UpdateTime, SwitchTime, Reps, Mean1, Cont1, Mean2, Cont2
```
* **`UpdateTime`** *(int, ms)*: Noise frame interval in milliseconds.
* **`SwitchTime`** *(int, ms)*: Duration in milliseconds before toggling to the alternative distribution.
* **`Reps`** *(int)*: Total number of distribution switch cycles.
* **`Mean1`, `Cont1`**: Mean ($0.0-1.0$) and Contrast/Std for State 1.
* **`Mean2`, `Cont2`**: Mean ($0.0-1.0$) and Contrast/Std for State 2.
* **Example:** `cs, 10, 3000, 10, 0.5, 0.05, 0.5, 0.25` *(Alternates every 3s between 5% and 25% contrast)*.

---

### 5. Frequency Chirp / Sweep (`fs`)
Performs an exponential frequency chirp starting from $F_{\text{min}}$, sweeping up to $F_{\text{max}}$, and returning to $F_{\text{min}}$.
```text
fs, Fmin, Fmax, SweepFactor, PhaseA, PhaseB, ContrastA, ContrastB
```
* **`Fmin` / `Fmax`** *(float, Hz)*: Minimum and maximum frequencies of the sweep.
* **`SweepFactor`** *(float)*: Rate of exponential frequency transition.
* **`PhaseA` / `PhaseB`** *(float, 0.0 – 1.0)*: Channel phase offsets.
* **`ContrastA` / `ContrastB`** *(float, 0.0 – 1.0)*: Channel contrast depths.

---

### 6. Timed Luminance Step (`sdt`)
Applies a static DC duty cycle for a specified duration before resetting to zero.
```text
sdt, DutyA, DutyB, Duration
```
* **`DutyA` / `DutyB`** *(float, 0.0 – 100.0)*: Output percentage for Channel A and Channel B.
* **`Duration`** *(int, ms)*: Duration in milliseconds.
* **Example:** `sdt, 50.0, 0.0, 1000` *(1-second 50% flash on Channel A)*.

---

### 7. Configuration & Utility Commands

#### Static Duty (`sd`)
```text
sd, DutyA, DutyB
```
Sets continuous DC output without timeout. Example: `sd, 0, 0` turns off both channels.

#### Gamma Calibration Sweep (`gc`)
```text
gc, Step, Wait, Reps
```
Iterates through luminance levels from $0\%$ to $100\%$ in increments of `Step` percent, holding each level for `Wait` milliseconds, repeated `Reps` times. Used with a power meter or photodiode.

#### Active Gamma LUT Selection (`agc`)
```text
agc, Index
```
Switches active gamma correction LUT bank (`1` = LUT 1, `2` = LUT 2).

#### Analog Calibration Stream (`ana`)
```text
ana
```
Streams continuous analog sensor values from pins A0 and A1 over Serial to monitor manual potentiometer adjustments. Send `done\n` to terminate stream.
