# Integration & Multi-Modal Synchronization

Complex experiments often require synchronized acquisition across multiple data streams: visual stimulation timestamps, high-speed eye tracking, inertial motion measurements, and neural electrophysiology or imaging.

!!! note "Work in Progress"
    This documentation section is currently a work in progress.

---

## Synchronization Architecture

```mermaid
flowchart TD
    subgraph Hardware["⚡ Hardware TTL Signals"]
        direction LR
        MCU_P4["<b>Microcontroller Pin 4</b><br/><i>Indicator / Frame Toggle</i>"]
        MCU_P5["<b>Microcontroller Pin 5</b><br/><i>Stimulus Active Gate</i>"]
    end

    subgraph Software["💻 Software Acquisition (LabStreamingLayer)"]
        direction LR
        Bonsai["<b>Bonsai Stimulus Engine</b><br/><i>LSL Outlet: Stimulus Events</i>"]
        Neon["<b>Pupil Labs Neon Tracker</b><br/><i>LSL Outlet: Gaze and Video</i>"]
        IMU["<b>Motion / IMU Stream</b><br/><i>LSL Outlet: 6-DOF Kinematics</i>"]
    end

    SyncBox["<b>Hardware DAQ / Sync Box</b><br/><i>Open Ephys / Intan / NI-DAQ</i>"]
    LabRecorder["<b>LabRecorder Multi-Stream Host</b><br/><i>Synchronized XDF Output File</i>"]

    MCU_P4 --> SyncBox
    MCU_P5 --> SyncBox
    Bonsai --> LabRecorder
    Neon --> LabRecorder
    IMU --> LabRecorder
    SyncBox -.->|Hardware Clocks| LabRecorder
```

---

## 1. Hardware TTL Synchronization Pins

The microcontroller outputs real-time TTL state flags to synchronize external acquisition systems (e.g. Open Ephys, Intan, National Instruments DAQ):
* **Pin 4 (Indicator Toggle):** Inverts logic state on every stimulus update, frame transition, or trial onset.
* **Pin 5 (Stimulus Status):** Drives `HIGH` throughout the active stimulus duration and returns to `LOW` during inter-trial intervals (ITIs).

---

## 2. Pupil Labs Neon Eye Tracking Integration

* **Workflow:** [`BonsaiCode/streamNeonDev.bonsai`](file:///d:/Code/LEDStimulation/BonsaiCode/streamNeonDev.bonsai)
* **Description:** Connects to the Pupil Labs Neon mobile eye-tracking stream over network API / RTSP / LSL to record:
  * 200 Hz binocular eye camera video and pupil diameter.
  * Real-time gaze coordinates in pixel and visual angle coordinates.
  * Scene camera video synchronized with optical stimulation onset.

---

## 3. Inertial Measurement Unit (IMU) Tracking

* **Firmware:** [`ArduinoCode/IMUscript_Task/`](file:///d:/Code/LEDStimulation/ArduinoCode/IMUscript_Task/)
* **Analysis:** [`HumanExpAnalysis/imu_analysis_dev.m`](file:///d:/Code/LEDStimulation/HumanExpAnalysis/imu_analysis_dev.m)
* **Application:** Records linear acceleration and angular velocity during mobile human or animal tasks to study vision during locomotion, head turns, or balance challenges.

---

## 4. LabRecorder & LabStreamingLayer (LSL)

* **Directory:** [`LabRecorderFiles/`](file:///d:/Code/LEDStimulation/LabRecorderFiles/)
* **Format:** Extensible Data Format (`.xdf`)
* **Benefit:** Ensures sub-millisecond software synchronization across heterogeneous devices (Bonsai, Pupil Labs Neon, IMUs, BioSemi / EEG amplifiers) via automatic clock-offset correction.
