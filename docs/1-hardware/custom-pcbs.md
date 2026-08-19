# Custom PCBs

This project includes custom Printed Circuit Boards (PCBs) designed in KiCad for driving multi-wavelength LED arrays and interfacing with the control electronics.

---

## 1. Standard LED Array Board (`LED_board`)

* **Location:** [`Hardware/LED_board/`](file:///d:/Code/LEDStimulation/Hardware/LED_board/)
* **Files:** `LED_board.kicad_sch`, `LED_board.kicad_pcb`, `LED_board.step`, `LED_board.stl`, `production/`
* **Status:** **Active**

### Description
The standard LED array board hosts multi-channel LED strings designed for uniform full-field illumination when mounted within the diffuser enclosure.
- **Topology:** Dual-channel parallel/series LED strings with dedicated current-limiting and connector interfaces.
- **Connectors:** Right-angled headers for modular wiring and daisy-chaining.
- **Mounting:** Screw holes aligned with the custom aluminium baseplate and MakerBeam framing.

---

## 2. Human Experiment LED Board (`LED_board_humans`)

* **Location:** [`Hardware/LED_board_humans/`](file:///d:/Code/LEDStimulation/Hardware/LED_board_humans/)
* **Status:** **Active / Application Specific**

### Description
Adapted layout customized for human psychophysics and visual presentation setups (e.g. specialized spacing or wearable diffuser mounting).

---

## 3. Constant Current Driver Board (`Constant_current_driver`)

* **Location:** [`Hardware/Constant_current_driver/`](file:///d:/Code/LEDStimulation/Hardware/Constant_current_driver/)
* **Files:** `Constant_current_driver.kicad_sch`, `Constant_current_driver.kicad_pcb`
* **Status:** ⚠️ **WIP / In Development (Not Active)**

> [!WARNING]
> This driver board is an experimental prototype and is **not currently active** in the standard experimental pipeline. Refer to the main [Circuit & Driver Design](file:///d:/Code/LEDStimulation/docs/1-hardware/circuit-and-driver.md) for current driving configurations.

---

## Fabrication & Production

Gerber files and drill maps for active boards are located in their respective `production/` subdirectories, ready for standard PCB fabrication (e.g. JLCPCB, PCBWay).
