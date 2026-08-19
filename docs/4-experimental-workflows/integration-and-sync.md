# Integration & Multi-Modal Synchronization

Complex experiments often require synchronized acquisition across multiple data streams: visual stimulation timestamps, high-speed eye tracking, inertial motion measurements, and neural electrophysiology or imaging.

---

## Synchronization Architecture

```mermaid
graph TD
    subgraph Hardware Signals
        MCU[Microcontroller Pin 4: Indicator Pin] --> SyncBox[Hardware Sync / DAQ]
        MCU[Microcontroller Pin 5: Active Status] --> SyncBox
    end

    subgraph Software Synchronization - LabStreamingLayer
        Bonsai[Bonsai Stimulus Orchestrator] -- LSL Outlet: Stimulus Events --> LSL_Router[LSL Network Router]
        Neon[Pupil Labs Neon Eye Tracker] -- LSL Outlet: Gaze & Pupil Video --> LSL_Router
        IMU[IMU / Motion Sensor Stream] -- LSL Outlet: Accelerometer/Gyro --> LSL_Router
        
        LSL_Router --> LabRecorder[LabRecorder (XDF Multi-Stream File)]
    end
    
    SyncBox -. Hardware Clocks .-> LabRecorder
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
