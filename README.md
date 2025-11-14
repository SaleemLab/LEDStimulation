# LED stimulator


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



# Stimuli

| Stimulus Category        | Stimulus Type                 | Chroma              | Luminance            | Parameters / Notes                                  |
| :----------------------- | :---------------------------- | :------------------ | :------------------- | :-------------------------------------------------- |
| **LED Stimuli** | Sine wave flicker             | green, UV, green+UV | $2 \times 10^3$, $2 \times 10^2$ |                                                     |
|                          | Sine wave envelope            | green, UV           | $2 \times 10^3$, $2 \times 10^2$ | (2p + verify on ephys)                             |
|                          | Full-field flash              | UV, green, UV+green | $2 \times 10^3$, $2 \times 10^2$ | (i.e. square wave flicker)                                |
|                          | White noise                   | green, UV, (maybe combined?) | $2 \times 10^3$     |                                                     |
|                          | Frozen white noise            | green, UV, (maybe combined?) | $2 \times 10^3$     |                                                     |
|                          | Contrast Adaptation           | green, UV                     | $2 \times 10^3$     | Periods: 4,8,16,32s    Contrast levels: 10 vs 100?                                                |
|                          | Luminance Adaptation            | green, UV,          | $2 \times 10^3$     |       $2 \times 10^1$ to  $2 \times 10^3$, 5% or 10% contrast super-imposed                                        |
