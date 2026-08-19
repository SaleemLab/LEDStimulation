# Bill of Materials (BOM)

This Bill of Materials is categorized by subsystem. Sourcing links and part references are provided below.

---

## 1. Full-Field Panel Enclosure & Mechanics

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Chassis** | MakerBeam XL 200mm (Black Anodized) | [Technobots Online](https://www.technobotsonline.com/makerbeamxl-200mm-long-black-anodised-beam-threaded.html) | 8 | Main cube frame |
| **Chassis** | MakerBeam XL 50mm (Black Anodized) | [Technobots Online](https://www.technobotsonline.com/makerbeamxl-50mm-long-black-anodised-beam-threaded.html) | 4 | Vertical riser beams |
| **Chassis** | MakerBeam XL Corner Cubes (Pack of 12) | [Technobots Online](https://www.technobotsonline.com/makerbeam-xl-black-corner-cubes-pack-of-12.html) | 8 cubes | Corner joints |
| **Fasteners** | M3 Button Head Socket Screws (6mm / 10mm) | [Technobots Online](https://www.technobotsonline.com/button-head-socket-cap-stainless-screw-m3x6mm-pk-100.html) | ~200 | M3 hardware |
| **Baseplate** | Custom Aluminium Baseplate | [`Hardware/Enclosure/Baseplate_for_PCBs_v3.stp`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/Baseplate_for_PCBs_v3.stp) | 1 | CNC milled aluminium |
| **Diffusion** | Corner Diffusion Base (Screws + Magnets) | [`Hardware/Enclosure/box_corner_screws_and_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_corner_screws_and_magnets.stl) | 4 | 3D printed |
| **Diffusion** | Edge Diffusion Base (Screws + Magnets) | [`Hardware/Enclosure/box_long_edge_screws_and_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_long_edge_screws_and_magnets.stl) | 4 | 3D printed |
| **Diffusion** | Corner Diffusion Magnetic Clamp | [`Hardware/Enclosure/box_corner_just_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_corner_just_magnets.stl) | 4 | 3D printed |
| **Diffusion** | Edge Diffusion Magnetic Clamp | [`Hardware/Enclosure/box_long_edge_just_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_long_edge_just_magnets.stl) | 4 | 3D printed |
| **Magnets** | 5mm Neodymium Disc Magnets | [Amazon UK](https://www.amazon.co.uk/Magnet-Expert%C2%AE-5mm-thick-Neodymium/dp/B00TACMMP0/) | 8+ | Embedded in 3D prints |
| **Diffuser** | Optical Diffuser Sheet (e.g. holographic/opal) | Lab supplier / Edmund Optics / Thorlabs | 1 | Sized to 200x200mm |

---

## 2. Optoelectronics & LED Boards

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **LEDs** | Green LEDs (e.g. 525nm peak) | [LCSC Electronics](https://www.lcsc.com/) | As per PCB | M-cone / Rod stimulation |
| **LEDs** | UV LEDs (e.g. 365–405nm peak) | [LCSC Electronics](https://www.lcsc.com/) | As per PCB | S-cone stimulation |
| **PCB** | Custom Dual-Channel LED PCB | [`Hardware/LED_board/`](file:///d:/Code/LEDStimulation/Hardware/LED_board/) | 1+ | Standard panel array |
| **PCB** | Human Miniature LED PCB | [`Hardware/LED_board_humans/`](file:///d:/Code/LEDStimulation/Hardware/LED_board_humans/) | 1 | Wearable glasses array |
| **Hardware** | M3 Nylon Spacers (2mm) | Standard supplier | 24 | PCB standoff isolation |
| **Hardware** | M3 x 10mm Flanged Button Screws & Nuts | [Accu](https://www.accu.co.uk/flanged-button-screws/8595-SSBF-M3-10-A2) | 24 | PCB mounting |
| **Headers** | Right-Angled Pin Headers (Male/Female) | Standard supplier | – | Modular inter-board wiring |

---

## 3. Control & Driver Electronics

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Microcontroller** | Arduino Leonardo (ATmega32u4) | Arduino Official / Authorized Distributor | 1 | Default controller |
| **Microcontroller** | Teensy 4.1 (ARM Cortex-M7) *(Optional)* | PJRC / SparkFun | 1 | High-speed DDS variant |
| **Driver Stage** | Power N-Channel MOSFETs | Standard supplier | 2+ | Logic-level or low-$R_{DS(on)}$ |
| **Protection** | Gate Resistors (100 $\Omega$) & Pull-downs (10 k$\Omega$) | Standard supplier | 2+ | Switching stability |
| **Protection** | Snubber Capacitors (100 nF) | Standard supplier | 2+ | High-frequency ringing filter |
| **Potentiometers** | Precision 10-turn / Trim Potentiometers | Standard supplier | 2 | Analog gate voltage scaling |

---

## 4. Wearable & Human Experiment Accessories

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Eye Tracker** | Pupil Labs Neon Eye Tracker | [Pupil Labs](https://pupil-labs.com/products/neon/) | 1 | High-speed mobile gaze tracking |
| **Mounting** | 3D Printed Glasses Boom Arm | [`Hardware/glasses/boom_arm.stl`](file:///d:/Code/LEDStimulation/Hardware/glasses/boom_arm.stl) | 1 | Mounts to Neon frame |
| **Diffuser Housing**| Miniature Glasses Diffuser Housing | [`Hardware/glasses/diffuser_fullipt.ipt`](file:///d:/Code/LEDStimulation/Hardware/glasses/diffuser_fullipt.ipt) | 1 | Sits in front of eye |
| **Motion Sensor** | IMU / Accelerometer *(Optional)* | Standard supplier | 1 | Body/gait synchronization |

---

## 5. Calibration Equipment

| Category | Part Description | Source / Reference | Notes |
| :--- | :--- | :--- | :--- |
| **Power Meter** | Optical Power Meter Console (e.g. PM100D) | Thorlabs | Absolute power measurements |
| **Photodiode** | Calibrated UV/Visible Photodiode Sensor | Thorlabs (e.g. S120VC / S121C) | Calibrated spectral response |
| **Spectrometer** | Spectroradiometer *(Optional)* | Ocean Optics / Thorlabs | Exact LED emission spectra |
