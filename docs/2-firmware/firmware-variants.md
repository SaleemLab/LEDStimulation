# Firmware Variants

The repository contains several firmware variants tailored for different microcontrollers, signal generation techniques, and experimental rigs.

---

## Variant Overview

| Sketch Directory | Target MCU | Core Engine | Primary Use Case | Status |
| :--- | :--- | :--- | :--- | :--- |
| [`StimulusControlViaSerial_Leonardo_v8`](file:///d:/Code/LEDStimulation/ArduinoCode/StimulusControlViaSerial_Leonardo_v8/) | Arduino Leonardo (ATmega32u4) | Timer1 PWM + Timer3 CTC ISR | General-purpose visual stimulator, full serial API | **Active / Canonical** |
| [`StimulusControlViaSerial_Leonardo_DDS_8bit_mouseRig`](file:///d:/Code/LEDStimulation/ArduinoCode/StimulusControlViaSerial_Leonardo_DDS_8bit_mouseRig/) | Arduino Leonardo (ATmega32u4) | Direct Digital Synthesis (DDS) phase accumulator | Mouse rig experiments requiring ultra-fast continuous phase updates | **Active / Rig-Specific** |
| [`StimulusControlViaSerial_Leonardo_DDS_8bit_2freq`](file:///d:/Code/LEDStimulation/ArduinoCode/StimulusControlViaSerial_Leonardo_DDS_8bit_2freq/) | Arduino Leonardo (ATmega32u4) | Dual-Frequency DDS | Simultaneous dual-frequency carrier modulation | **Experimental** |
| [`LEDStimController_Teensy41_DDS`](file:///d:/Code/LEDStimulation/ArduinoCode/LEDStimController_Teensy41_DDS/) | Teensy 4.1 (ARM Cortex-M7 @ 600MHz) | High-speed DDS + Hardware Timers | Ultra-high sample rate and multi-channel expansion | **Active / Advanced** |
| [`IMUscript_Task`](file:///d:/Code/LEDStimulation/ArduinoCode/IMUscript_Task/) | Arduino Leonardo / Compatible | Serial IMU Data Acquisition | Movement / gait synchronization in human paradigms | **Active / Auxiliary** |

---

## 1. Leonardo v8 (`StimulusControlViaSerial_Leonardo_v8`)
* **Target:** Arduino Leonardo (ATmega32u4)
* **Architecture:** 16-bit Phase Correct PWM (Timer1) with dynamic interrupt callback pacing (Timer3).
* **Highlights:** Full implementation of all serial stimulus types (`s`, `wn`, `fwn`, `cs`, `fs`, `sdt`), LUT gamma mapping, and potentiometer analog reading. Recommended for standard lab setups.

---

## 2. Leonardo DDS (`StimulusControlViaSerial_Leonardo_DDS_8bit_mouseRig`)
* **Target:** Arduino Leonardo (ATmega32u4)
* **Architecture:** Direct Digital Synthesis (DDS) using fixed-point phase accumulators for sine generation.
* **Highlights:** Allows instant, glitch-free frequency and phase transitions during live experimental workflows.

---

## 3. Teensy 4.1 DDS (`LEDStimController_Teensy41_DDS`)
* **Target:** Teensy 4.1 (600 MHz ARM Cortex-M7)
* **Architecture:** Hardware FlexPWM / IntervalTimers combined with 32-bit phase accumulators.
* **Highlights:** Significantly higher processing headroom, nanosecond timing precision, and capability for high-rate continuous modulation without MCU bottleneck.
