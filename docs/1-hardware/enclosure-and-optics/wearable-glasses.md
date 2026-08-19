# Wearable Glasses Setup (Pupil Labs Neon Integration)

For human visual psychophysics, mobile experiments, and simultaneous eye-tracking paradigms, a miniature stimulation and diffusion assembly is mounted directly to the **Pupil Labs Neon** eye-tracking frame.

---

## Hardware & 3D Models

* **CAD Location:** [`Hardware/glasses/`](file:///d:/Code/LEDStimulation/Hardware/glasses/)
* **Key Files:**
  * `boom_arm.stl` / `boom_arm.ipt` — Custom 3D printed boom arm that attaches securely to the Neon eye-tracker frame.
  * `diffuser_fullipt.ipt` — Miniature optical diffuser housing that holds the stimulus LED array in front of the subject's visual field without obstructing the eye camera's gaze-tracking optics.
  * `Neon_JAN_frame_reference/` — Reference frame geometries for fit and alignment validation.

---

## Assembly & Alignment

```mermaid
graph LR
    NeonFrame[Pupil Labs Neon Frame] --> BoomArm[Custom 3D Printed Boom Arm]
    BoomArm --> DiffuserHousing[Miniature Diffuser Housing]
    DiffuserHousing --> HumanPCB[Human LED Board]
    DiffuserHousing --> Eye[Subject Visual Field]
    NeonFrame -.-> GazeTracking[Unobstructed Eye-Tracking Cameras]
```

1. **Boom Arm Attachment:** The boom arm securely clips/fastens to the Neon eye-tracker frame without altering gaze-calibration geometry.
2. **Diffuser & LED Insertion:** The miniature human LED board (`LED_board_humans`) seats into the diffuser housing, placing an evenly diffused light source in the peripheral or full visual field.
3. **Cable Routing:** Lightweight, flexible silicone wiring runs along the glasses temple to prevent torque or movement artifacts during natural head motion and walking tasks.
