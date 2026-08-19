# Bill of Materials (BOM)

This Bill of Materials is categorized by subsystem. Sourcing links, specific part numbers, and references are provided below.

---

## 1. Full-Field Panel Enclosure & Mechanics

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Chassis** | MakerBeam XL 200mm (Black Anodized) | [Technobots Online](https://www.technobotsonline.com/makerbeamxl-200mm-long-black-anodised-beam-threaded.html) | 8 | Main cube frame |
| **Chassis** | MakerBeam XL 50mm (Black Anodized) | [Technobots Online](https://www.technobotsonline.com/makerbeamxl-50mm-long-black-anodised-beam-threaded.html) | 4 | Vertical riser beams |
| **Chassis** | MakerBeam XL Corner Cubes (Pack of 12) | [Technobots Online](https://www.technobotsonline.com/makerbeam-xl-black-corner-cubes-pack-of-12.html) | 8 cubes | Corner joints |
| **Fasteners** | M3 Button Head Socket Screws (6mm / 10mm) | [Technobots Online](https://www.technobotsonline.com/button-head-socket-cap-stainless-screw-m3x6mm-pk-100.html) | ~200 | Frame assembly & baseplate mounting |
| **Fasteners** | MakerBeam XL M3 T-Slot / Square Nuts | [Technobots Online](https://www.technobotsonline.com/) | 50+ | Slides into 15x15mm extrusion slots |
| **Baseplate** | Custom Aluminium Baseplate | [`Hardware/Enclosure/Baseplate_for_PCBs_v3.stp`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/Baseplate_for_PCBs_v3.stp) | 1 | CNC milled / laser cut aluminium heatsink plate |
| **Diffusion Base** | Corner Diffusion Base (Screws + Magnets) | [`Hardware/Enclosure/box_corner_screws_and_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_corner_screws_and_magnets.stl) | 4 | 3D printed |
| **Diffusion Base** | Edge Diffusion Base (Screws + Magnets) | [`Hardware/Enclosure/box_long_edge_screws_and_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_long_edge_screws_and_magnets.stl) | 4 | 3D printed |
| **Diffusion Clamp** | Corner Diffusion Magnetic Clamp | [`Hardware/Enclosure/box_corner_just_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_corner_just_magnets.stl) | 4 | 3D printed |
| **Diffusion Clamp** | Edge Diffusion Magnetic Clamp | [`Hardware/Enclosure/box_long_edge_just_magnets.stl`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/box_long_edge_just_magnets.stl) | 4 | 3D printed |
| **Magnets** | 5mm $\times$ 2–3mm Neodymium Disc Magnets | [Amazon UK](https://www.amazon.co.uk/Magnet-Expert%C2%AE-5mm-thick-Neodymium/dp/B00TACMMP0/) | 16 | 8 clamp positions $\times$ 2 magnets (base + clamp) |
| **Diffuser** | Tracing paper works well for diffusing and passing through the UV / visible light spectrum. Larger panels may benefit from multiple sheets. | 1 | Sized to 200 $\times$ 200 mm |
| **3D Print Set** | Complete Enclosure 3D Print Package | [`Hardware/Enclosure/UM3_all3dprintedparts.3mf`](file:///d:/Code/LEDStimulation/Hardware/Enclosure/UM3_all3dprintedparts.3mf) | 1 set | Combined Ultimaker / Bambu 3MF package |

---

## 2. Optoelectronics & LED Boards

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Green LEDs** | Xinglight XL-2835UGC-02 (525nm peak, SMD 2835) | [LCSC (C2843870)](https://www.lcsc.com/product-detail/Light-Emitting-Diodes-LED_XINGLIGHT-XL-2835UGC-02_C2843870.html) | 6 per board | M-cone / Rod stimulation (2 series banks of 3) |
| **UV LEDs** | Xinglight XL-2835UVA-02 (395–405nm peak, SMD 2835) | [LCSC (C2843877)](https://www.lcsc.com/product-detail/Light-Emitting-Diodes-LED_XINGLIGHT-XL-2835UVA-02_C2843877.html) | 6 per board | S-cone stimulation (2 series banks of 3) |
| **PCB** | Custom Dual-Channel LED PCB | [`Hardware/LED_board/`](file:///d:/Code/LEDStimulation/Hardware/LED_board/) | 1 to 4 | Standard panel array (JLCPCB production files in `production/`) |
| **PCB (3-Colour)**| 3-Colour LED PCB (`LED_board_humans`) | [`Hardware/LED_board_humans/`](file:///d:/Code/LEDStimulation/Hardware/LED_board_humans/) | – | *WIP* – Experimental 3-colour panel PCB |
| **Connectors** | 2-pin 2.54mm Right-Angled Female Sockets & Male Headers | Standard supplier | 4 per board | `PinSocket_1x02_P2.54mm_Horizontal` for LED bank wiring |
| **Hardware** | M3 Nylon Spacers (2mm) | Standard supplier | 4–16 | PCB standoff thermal and electrical isolation |
| **Hardware** | M3 $\times$ 10mm Flanged Button Screws & Nuts | [Accu](https://www.accu.co.uk/flanged-button-screws/8595-SSBF-M3-10-A2) | 4–16 | PCB mounting to baseplate |

---

## 3. Control & Driver Electronics

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Microcontroller** | Arduino Leonardo (ATmega32u4, 5V / 16MHz) | Arduino Official / Authorized Distributor | 1 | Default controller (Timer 1/3 high-speed PWM) |
| **Microcontroller** | Teensy 4.1 (ARM Cortex-M7, 600MHz) *(Optional)* | PJRC / SparkFun | 1 | High-speed DDS variant |
| **Transistor Switch** | NPN Bipolar Junction Transistor (e.g. onsemi 2N3904, TO-92) | Standard supplier (Mouser / Farnell / Digikey) | 1 per channel | Low-side switch ($40\,\text{V}$, $200\,\text{mA}$ max) |
| **Potentiometers** | Multi-Turn Cermet Trimpot (e.g. 500 $\Omega$ – 10 k$\Omega$, Bourns 3296W) | Standard supplier | 1 per channel | Master peak current / luminance control (cermet is critical for PWM stability) |
| **Fixed Resistors** | $R_{\text{fixed}}$ Series Current-Limiting Resistors (51 $\Omega$ – 330 $\Omega$, 1/4W) | Standard supplier | 1 per LED line | Protects LEDs and sets current ceiling when trimpot is at $0\,\Omega$ |
| **Base Resistors** | $R_{\text{base}}$ Base Resistors (100 $\Omega$, 1/4W) | Standard supplier | 1 per channel | Limits base switching current from Arduino pin |
| **Pull-down Resistors** | $R_{\text{pull-down}}$ Pull-Down Resistors (10 k$\Omega$, 1/4W) | Standard supplier | 1 per channel | Keeps transistor OFF during boot / floating GPIO |
| **Power Supply** | $+12\,\text{V}$ DC Power Supply (1A–2A, Regulated) | Standard supplier | 1 | External power rail for multi-LED series strings |
| **Prototyping** | Solderless Breadboard or Stripboard / Perfboard | Standard supplier | 1 | Benchtop circuit assembly |
| **Wiring** | USB-A to Micro-USB Cable & Jumper Wires | Standard supplier | As needed | Serial communication, power & interconnections |

---

## 4. Wearable & Human Experiment Accessories

| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Eye Tracker** | Pupil Labs Neon Eye Tracker | [Pupil Labs](https://pupil-labs.com/products/neon/) | 1 | High-speed mobile gaze tracking system |
| **Mounting** | 3D Printed Glasses Boom Arm | [`Hardware/glasses/boom_arm.stl`](file:///d:/Code/LEDStimulation/Hardware/glasses/boom_arm.stl) | 1 | Custom clip mounting to Neon frame |
| **Diffuser Housing**| Miniature Glasses Diffuser Housing | [`Hardware/glasses/diffuser_fullipt.ipt`](file:///d:/Code/LEDStimulation/Hardware/glasses/diffuser_fullipt.ipt) | 1 | Positions diffuser in visual field |
| **Optics** | Miniature Stimulus LED (UV / Green) | Standard supplier / discrete LEDs | 1–2 | Inserted into glasses diffuser housing |
| **Wiring** | Ultra-flexible Silicone Hookup Wire (30 AWG) | Standard supplier | 1–2 m | Lightweight routing along glasses temple |
| **Motion Sensor** | IMU / Accelerometer *(Optional)* | Standard supplier | 1 | Body/gait synchronization |

---

## 5. Calibration Equipment

!!! important "Only ONE Calibration Setup Required"
    **Only ONE calibrated optical measurement setup is required** (Options A, B, or C below). You do **not** need all of them. Choose the option best suited to your available lab equipment and acquisition requirements.

### Option A: Standard Optical Power Meter (Recommended)
| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Power Meter** | Optical Power Meter Console (e.g. Thorlabs PM100D, PM100USB, or PM400) | [Thorlabs](https://www.thorlabs.com/) | 1 | Measures absolute optical power ($\mu\text{W}$ / $\text{nW}$) |
| **Sensor** | Calibrated UV/Visible Photodiode Sensor (e.g. Thorlabs S120VC, S121C, S130VC, or S120C) | [Thorlabs](https://www.thorlabs.com/) | 1 | Wavelength-calibrated response curve ($200\text{--}1100\,\text{nm}$) |

### Option B: High-Speed Amplified Photodiode (For Dynamic / Temporal Profiling)
| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Photodiode** | Amplified Si Photodetector (e.g. Thorlabs PDA100A2) | [Thorlabs](https://www.thorlabs.com/) | 1 | High bandwidth photodiode for fast temporal/PWM waveform verification |
| **Acquisition** | Oscilloscope or National Instruments / LabJack DAQ | Standard lab equipment | 1 | Digitizes amplified analog photodiode voltage output |

### Option C: Calibrated Spectroradiometer / Spectrometer
| Category | Part Description | Source / Reference | Quantity | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **Spectrometer** | Calibrated Fiber Spectroradiometer (e.g. Ocean Optics USB2000+ / STS, Thorlabs CCS200) | Ocean Optics / Thorlabs | 1 | Directly measures full spectral irradiance distributions ($\mu\text{W/cm}^2\text{/nm}$) |
