# Full-Field LED Visual Stimulator

This project implements a high-precision, large-field LED stimulator designed for visual neuroscience applications. It drives dual-channel LED arrays (e.g., Green for M-cones/rods and UV for S-cones) to produce precise temporal stimulation patterns such as sinusoidal flickers, Gaussian white noise (including contrast and mean switching adaptation paradigms), luminance steps and frequency chirps.

This ReadME is a WIP.

## System Overview

The stimulator runs on an **Arduino Leonardo (ATmega32u4)**. It utilizes a hybrid control mechanism to manage LED brightness:
1.  **Fast PWM Signal:** Generates temporal patterns using high-frequency digital switching.
2.  **Analog Dimming:** Uses a potentiometer to scale the gate voltage, forcing the driving MOSFETs to operate in the linear (ohmic) region. This acts as an adjustable current limiter for setting peak luminance.

## Hardware Setup

### Circuit Design
The system uses a low-side switching topology.
*   **LED Load:** Parallel strings of 3 LEDs, each with a current-limiting resistor.
*   **Switching:** N-channel MOSFETs (one per chromatic channel).
*   **Protection:**
    *   **Snubber Capacitor:** Placed across the MOSFET (Source-Drain) to dampen high-frequency ringing.
    *   **Gate Resistor:** 100$\Omega$ between Arduino and Gate to limit inrush current.
    *   **Pull-down Resistor:** 10k$\Omega$ from Gate to Ground to ensure the MOSFET stays off during boot.
 
[ciruit schematic]

### Pinout Configuration
Based on the firmware configuration, connect your hardware to the Arduino Leonardo as follows:

| Component | Arduino Pin | Description |
| :--- | :--- | :--- |
| **Channel A (LEDs)** | **Pin 9** | PWM Output (Driven by Timer1) |
| **Channel B (LEDs)** | **Pin 10** | PWM Output (Driven by Timer1) |
| **Indicator** | **Pin 4** | Toggles state to indicate stimulus change events |
| **Stimulus Status** | **Pin 5** | HIGH when stimulus is active, LOW otherwise |
| **Sensors** | **A0 / A1** | For setting the potentiometer values |

## Firmware Features

*   **Dual-Channel Control:** Independent modulation of Channel A and Channel B.
*   **High-Resolution PWM:** Uses 16-bit Phase Correct PWM (Timer1) to prevent visible digital stepping.
*   **Precise Timing:** Waveform updates are driven by hardware interrupts (Timer3) rather than software `delay()` loops, ensuring stable temporal frequencies.
*   **Gamma Correction:** Implements Look-Up Tables (LUTs) stored in PROGMEM to linearize LED output.
*   **Fast RNG:** Uses `xorshift32` for high-speed random number generation required by white noise stimuli, approximating a Gaussian distribution via the Central Limit Theorem.

## Installation

1.  Open the firmware source code in the Arduino IDE.
2.  Select **Arduino Leonardo** as the board.
3.  Upload the sketch.
4.  Open the Serial Monitor or your preferred interface software (e.g., Bonsai).

## Practical User Guide: Serial API

The device is controlled via a serial interface. Send commands as comma-separated ASCII strings terminated by a carriage return (`\r`) or newline.

**Baud Rate:** `115200`.

### Stimulus Commands

#### 1. Sinusoidal Flicker (`s`)
Generates one or two channel sinusoidal flicker. Stimulus duration, frequency, phase and contrast can all be set.
```text
s, Duration, Freq, PhaseA, PhaseB, ContrastA, ContrastB
```
*   **Duration:** Stimulus length in milliseconds.
*   **Freq:** Temporal frequency in Hz.
*   **PhaseA/B:** Phase offset (0.0–1.0, where 1.0 = $2\pi$).
*   **ContrastA/B:** Modulation depth (0.0–1.0).

#### 2. Gaussian White Noise (`wn`)
Generates random luminance values based on a Gaussian distribution.
```text
wn, Duration, UpdateTime, MeanFrac, StdFrac
```
*   **UpdateTime:** Duration (ms) of each noise "frame".
*   **MeanFrac:** Center luminance level (0.0–1.0).
*   **StdFrac:** Standard deviation (sets the contrast).

#### 3. Frozen White Noise (`fwn`)
Plays a deterministic, repeatable sequence of Gaussian noise.
```text
fwn, Duration, UpdateTime, Reps, Seed
```
*   **Reps:** Number of times to loop the sequence.
*   **Seed:** Integer seed for the pseudo-random number generator.

#### 4. Switching White Noise (`cs`)
Alternates between two gaussian distributions (e.g., high vs. low variance, high vs low mean).
```text
cs, UpdateTime, SwitchTime, Reps, Mean1, Cont1, Mean2, Cont2
```
*   **SwitchTime:** Duration (ms) to stay in one state before switching.
*   **Mean/Cont:** Mean and Contrast parameters for State 1 and State 2.
*   **Reps:** number of distribution switches to perform.

#### 5. Frequency Chirp / Sweep (`fs`)
Performs an exponential frequency sweep (from Fmin to Fmax to Fmin)
```text
fs, Fmin, Fmax, SweepFactor, PhaseA, PhaseB, ContrastA, ContrastB
```
*   **Fmin/Fmax:** Start and End frequencies in Hz.
*   **SweepFactor:** Rate of change for the exponential sweep.

### 6. Timed luminance step
Performs a luminance step for a fixed duration before returning to baseline
```text
sdt, sdt, DutyA, DutyB, Duration
```
*   **DutyA/B:** Duty cycle value for each chromatic channel.


### Configuration Commands

| Command | Syntax | Description |
| :--- | :--- | :--- |
| **Set Static Duty** | `sd, DutyA, DutyB` | Sets DC output (0–100%) |
| **Gamma Calibration** | `gc, Step, Wait, Reps` | Steps through duty cycles to measure gamma response curve. |
| **Select LUT** | `agc, Index` | Switches Gamma LUTs. `1` = Chx1LUT, `2` = Chx2LUT. Can be autotoggled based on potentiometer settings|
| **Toggle Channel** | `useChA, 1/0` | Enables (`1`) or disables (`0`) Channel A/B output. |
| **Read Analog** | `ana` | Streams data from A0/A1 for configuring potentiometer values. Send "done" to stop. |

## Firmware Architecture & Timer Configuration

To achieve high-precision timing and smooth luminance gradients, the firmware directly manipulates the ATmega32u4 hardware timers.

### Timer1: High-Resolution PWM Generation
Timer1 is configured for **16-bit Phase Correct PWM** to drive the LED channels without visible digital stepping.

*   **Mode Configuration:** The Waveform Generation Mode (WGM) is set to Mode 10 (Phase Correct PWM, TOP = ICR1). This is achieved by setting the `WGM13` bit in `TCCR1B` and clearing `WGM12`, `WGM11`, and `WGM10`.
*   **Frequency Control:** The base frequency is determined by the `ICR1` register (Input Capture Register), which acts as the TOP value. The firmware dynamically calculates the optimal Prescaler (1, 8, 64, 256, or 1024) and `ICR1` value to approximate the target frequency (default ~7.68 kHz) while maximizing 16-bit resolution.
*   **Output Channels:**
    *   **Channel A (Pin 9):** Controlled by `OCR1A` (Output Compare Register). `COM1A1` is set to enable PWM output.
    *   **Channel B (Pin 10):** Controlled by `OCR1B`. `COM1B1` is set to enable PWM output.

### Timer3: Precision Stimulus Pacing
Timer3 drives the temporal updates of the stimulus (e.g., moving to the next point in a sine wave or generating the next white noise frame). It operates independently of the main loop.

*   **Mode Configuration:** Configured in **CTC Mode** (Clear Timer on Compare). The `WGM32` bit is set in `TCCR3B`.
*   **Interrupt-Driven:** The `OCIE3A` bit in `TIMSK3` is set to enable interrupts. When the timer counter matches the value in `OCR3A`, the `TIMER3_COMPA_vect` ISR fires.
*   **Dynamic Callbacks:** To support different stimulus types without rewriting the ISR, the firmware uses a function pointer (`timer3Callback`). The ISR executes whichever function (e.g., `sinewaveInterrupt`, `whiteNoiseInterrupt`) is currently assigned to this pointer.

### Summary of Register Settings

| Timer | Mode | Key Registers | Function |
| :--- | :--- | :--- | :--- |
| **Timer1** | 16-bit Phase Correct PWM | `TCCR1A`, `TCCR1B`, `ICR1` | Generates the carrier frequency for LED brightness. |
| **Timer3** | CTC (Clear Timer on Compare) | `TCCR3B`, `OCR3A`, `TIMSK3` | Triggers hardware interrupts to update stimulus waveforms. |



## Calibration
1) Determining potentiomter setting for specific luminance levels (can be measured in estimated photoisomerisations/s)
2) Generating gamma correction LUTs.

Thorlabs powermeter or photodiode. Must have UV sensitivity and calibration data.

## Synchronisation
Typical photodiode setup not possible since we are using a high frequency PWM which is ~16bit.


## BOM

| Category     | Part                                  | Link                                                                                               | Quantity  | Total Price (£) |
|--------------|----------------------------------------|---------------------------------------------------------------------------------------------------|----------|-----------------|
| **Enclosure**|                                        |                                                                                                   |          |                 |
|              | Custom aluminium baseplate            | [repo link](Hardware/Enclosure/Baseplate_for_PCBs_mirroredHole.stp)                                | 1        |  –               |
|              | 5mm aluminium side pieces             | –                                                                                                  | –        | –               |
|              | 200mm MakerBeam XL                    | [Link](https://www.technobotsonline.com/makerbeamxl-200mm-long-black-anodised-beam-threaded.html)  | 8       |             |
|              | 50mm MakerBeam XL                     | [Link](https://www.technobotsonline.com/makerbeamxl-50mm-long-black-anodised-beam-threaded.html)   | 4        |             |
|              | MakerBeam corner cubes (pack of 12)   | [Link](https://www.technobotsonline.com/makerbeam-xl-black-corner-cubes-pack-of-12.html)           | 8       |               |
|              | M3 button head socket screws (x100)   | [Link](https://www.technobotsonline.com/button-head-socket-cap-stainless-screw-m3x6mm-pk-100.html) | 200      |             |
|              | Diffusion base - corner               | [repo link](Hardware/Enclosure/box_corner_screws_and_magnets.stl)                                  | 4       | –               |
|              | Diffusion base - edge                 | [repo link](Hardware/Enclosure/box_long_edge_screws_and_magnets.stl)                               | 4       | –               |
|              | Diffusin clamp - corner               | [repo link](Hardware/Enclosure/box_corner_just_magnets.stl)                                        | 4        | –               |
|              | 5mm magnets                           | [Link](https://www.amazon.co.uk/Magnet-Expert%C2%AE-5mm-thick-Neodymium/dp/B00TACMMP0/ref=sr_1_6?s=kitchen)                      | 8        | –               |
|              |                                       | –                                                                                                  | –        | –               |
| **LEDs and connectors**|                             | -                                                                                                  |         |                 |
|              | Green LEDs                            | [LCSC link](https://www.lcsc.com/product-detail/C2843870.html?s_z=n_C2843870)              | –       | –               |
|              | UV LEDs                               | [LCSC link](https://www.lcsc.com/product-detail/C2843870.html?s_z=n_C2843870)              | –       | –               |
|              | Right-angled female pin headers       | –                                                                                                  | –        | –               |
|              | Right-angled male pin headers         | –                                                                                                  | –       | –               |
|              | red-black wire                        | –                                                                                                  | –      | –               |
|              | M3 nylon spacers, 2mm                 | –                                                                                                  | 24        | –               |
|              | M3 x 10mm Socket Flanged Button Screws | [Link](https://www.accu.co.uk/flanged-button-screws/8595-SSBF-M3-10-A2)                          | 24       | –               |
|              | M3 nuts                               | –                                                                                                  | 24       | –               |
|              |                                       | –                                                                                                  | –         | –               |
| **Control circuit**|                                 |                                                                                                    |         |                 |
|              | N-mosfet                              | –                                                                                                 | 2         | –               |
|              | Mosfet gate Resistors (100ohm)        | –                                                                                                  | 2        | –               |
|              | Mosfet source drain capacitor (100nF) | –                                                                                                  | 2        | –               |



