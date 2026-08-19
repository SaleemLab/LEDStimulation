# Circuit & Driver Design

> [!NOTE]
> **Status: Under Revision / In Transition**
> The hardware driving circuitry and control topology are currently undergoing active redesign. Detailed schematics, pinout assignments, gate-drive topologies, and component values will be populated here once the new revisions are finalized.

---

## Overview & Scope

The driver subsystem is responsible for:
- Delivering high-frequency modulated current to multi-channel LED arrays (e.g., UV, Green, or custom spectral configurations).
- Ensuring nanosecond-to-microsecond switching fidelity for high-frequency PWM and continuous waveform generation.
- Providing protective circuitry (gate protection, pulldowns, transient snubbing).
- Interfacing directly with microcontroller digital PWM timer outputs and hardware status lines.

---

## Planned Sections

- **Power Distribution & Supply**: Voltage rails, current budgeting, and power isolation.
- **Switching Topology**: Driver stage, gate control, and slew rate optimization.
- **Pinout Map**: Detailed mapping between MCU hardware pins (PWM channels, indicator toggles, status flags, and sensor lines) and the driver stage.
- **Wiring & Interconnect**: Cabling, connectors, and assembly safety.
