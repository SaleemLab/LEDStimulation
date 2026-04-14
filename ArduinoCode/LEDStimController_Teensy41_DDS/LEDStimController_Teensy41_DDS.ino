// ==========================================================================
//  LEDStimController_Teensy41_DDS.ino
//
//  Teensy 4.1 DDS-based LED Stimulation Controller
//  Ported from Arduino Leonardo (ATmega32U4) version
//
//  PWM:  12-bit @ 36,621 Hz on pins 2 (ChA) and 3 (ChB)
//  DDS:  32-bit phase accumulator, 4096-entry sine LUT
//  ISR:  IntervalTimer (PIT) with NVIC priority 16
//
//  Serial Commands (carriage return terminated):
//    s,dur,fA,fB,phA,phB,conA,conB     - Sinewave flicker
//    se,dur,fA,fB,envF,conA,conB        - Sine with contrast envelope
//    fs,fmin,fmax,sweep,phA,phB,cA,cB   - Continuous frequency sweep
//    sfs,f0,f1,step,cyc,phA,phB,cA,cB   - Stepped frequency sweep
//    sd,dutyA,dutyB                      - Set static duty cycle (%)
//    sdt,dutyA,dutyB,durationMs          - Set duty cycle for duration
//    gc,step,wait,nReps                  - Cycle through duty cycles
//    useChA,0/1                          - Enable/disable channel A
//    useChB,0/1                          - Enable/disable channel B
//    agc,index                           - Select gamma LUT (1 or 2)
//    ana                                 - Read analog calibration voltages
//    bright,A/B,percent                  - Set brightness scale (0-100%)
//    loadlut,A/B,v0,v1,...,v4095         - Upload empirical gamma LUT
// ==========================================================================

#include <math.h>
#include <IntervalTimer.h>

// ========================== CONFIGURATION ==========================

#define PWM_PIN_A   2       // FlexPWM4, submodule 2, channel A
#define PWM_PIN_B   3       // FlexPWM4, submodule 2, channel B
#define PWM_FREQ    36621.09f  // 150 MHz / 4096 — exact 12-bit resolution
#define PWM_BITS    12
#define PWM_MAX     4095    // (1 << PWM_BITS) - 1

#define PIN_CYCLE_INDICATOR  4   // Toggles on carrier cycle completion
#define PIN_STIM_ON          5   // HIGH during stimulation

#define TABLE_SIZE  4096    // Sine wavetable entries (matches PWM resolution)
#define TARGET_RECONFIG_INTERVAL_US 1000L  // Frequency sweep update interval

// --- Fade & Envelope ---
#define FADE_DURATION_MS    50.0f    // Fade-in/out time in milliseconds (0 = disabled)
#define FADE_LUT_SIZE       64
#define FADE_LUT_MAX        (FADE_LUT_SIZE - 1)

// --- Fixed-point scale for 12-bit ---
#define FIXED_SCALE     4096
#define FIXED_SHIFT     12

// ========================== GAMMA LUTs ==========================
// Two hardcoded LUT slots — paste empirical calibration data here,
// or upload at runtime via serial "loadlut" command.

// LUT 1: Linear (identity mapping)
uint16_t gammaLUT1[PWM_MAX + 1];

// LUT 2: Placeholder for empirical data (initialised to linear)
uint16_t gammaLUT2[PWM_MAX + 1];

// Active gamma LUT pointers — volatile because the ISR reads them and the main loop writes them
uint16_t * volatile currentGammaLUT_A;
uint16_t * volatile currentGammaLUT_B;

// Post-LUT brightness scaling (0–4096 = 0%–100%)
volatile uint16_t brightnessScaleA = FIXED_SCALE;
volatile uint16_t brightnessScaleB = FIXED_SCALE;

// Last pre-gamma value written to each channel — used to immediately re-apply
// the output when brightness scale changes (e.g. from the bright,A,50 command)
volatile uint16_t lastPwmValA = 0;
volatile uint16_t lastPwmValB = 0;

// ========================== PWM VARIABLES ==========================

uint16_t maxBrightness;   // Artificial TOP for flexible brightness control
uint16_t midBrightness;   // maxBrightness / 2

float actualPWMFreq = PWM_FREQ;  // Hardware-accurate frequency for DDS calculations

// ========================== CHANNEL SELECTION ==========================

// volatile: written by main loop (serial command), read by ISR
volatile bool useChA = true;
volatile bool useChB = true;

// ========================== SERIAL ==========================

// numChars covers all normal commands. The longest (sfs, 9 params) is well under 50 chars.
// For loadlut, always use stream mode (send just "loadlut,A" then values separately);
// inline mode is limited to ~5-6 test values due to this 50-char buffer.
const byte numChars = 50;
char receivedChars[numChars];
bool newData = false;

// ========================== SINE WAVETABLE ==========================

uint16_t sineWaveTable[TABLE_SIZE];

// ========================== DDS VARIABLES ==========================

volatile uint32_t phaseAccumulatorA = 0;
volatile uint32_t phaseAccumulatorB = 0;
volatile uint32_t phaseIncrementA = 0;
volatile uint32_t phaseIncrementB = 0;

volatile uint32_t envAccumulator = 0;
volatile uint32_t envIncrement = 0;

volatile uint16_t contrastMultIntA = FIXED_SCALE;
volatile uint16_t contrastMultIntB = FIXED_SCALE;

// Tracks full sine wave cycles (for frequency sweeps)
volatile uint32_t completedCycles = 0;

// ========================== RAISED COSINE ENVELOPE ==========================

uint16_t raisedCosineLUT[FADE_LUT_SIZE];
volatile unsigned long currentTick = 0;         // Master clock (ISR tick counter)
volatile unsigned long fadeInterrupts = 0;
volatile unsigned long fadeOutStartInterrupt = 0;
volatile uint32_t envStep = 0;

// ========================== INTERVAL TIMER ==========================

IntervalTimer ddsTimer;

// ========================== DDS PHASE INCREMENT ==========================

// Pre-computed scale factor: 2^32 / actualPWMFreq
// 2^32 = 4294967296 — the full range of the 32-bit phase accumulator.
// One full sine cycle (0 to 2π) maps onto the full 0 to 2^32 range.
// When the accumulator overflows past 2^32, it wraps — giving continuous oscillation.
static float phaseScale;

uint32_t calcPhaseInc(float freq) {
  return (uint32_t)(freq * phaseScale);
}

// ========================== PWM OUTPUT ==========================

// Apply gamma correction + post-LUT brightness scaling, then write to PWM.
// NOTE: analogWrite() is called from ISR context (sinewaveInterrupt etc.).
// On Teensy 4.1 this works reliably at 36.6 kHz — the function is fast and
// the FlexPWM double-buffer ensures glitch-free updates. If issues arise,
// replace with direct IMXRT_FLEXPWM4.SM[2].VAL3 / VAL5 register writes.
inline void setChA(uint16_t value) {
  lastPwmValA = value;  // Track for immediate re-apply on brightness change
  uint16_t gammaVal = currentGammaLUT_A[value];
  uint16_t scaled = ((uint32_t)gammaVal * brightnessScaleA) >> FIXED_SHIFT;
  analogWrite(PWM_PIN_A, scaled);
}

inline void setChB(uint16_t value) {
  lastPwmValB = value;  // Track for immediate re-apply on brightness change
  uint16_t gammaVal = currentGammaLUT_B[value];
  uint16_t scaled = ((uint32_t)gammaVal * brightnessScaleB) >> FIXED_SHIFT;
  analogWrite(PWM_PIN_B, scaled);
}

// ========================== SETUP ==========================

void setup() {
  // --- PWM Pins ---
  pinMode(PWM_PIN_A, OUTPUT);
  pinMode(PWM_PIN_B, OUTPUT);
  analogWriteFrequency(PWM_PIN_A, PWM_FREQ);  // Sets both pins (same FlexPWM submodule)
  analogWriteResolution(PWM_BITS);

  // --- Indicator Pins ---
  pinMode(PIN_CYCLE_INDICATOR, OUTPUT);
  pinMode(PIN_STIM_ON, OUTPUT);
  digitalWriteFast(PIN_CYCLE_INDICATOR, LOW);
  digitalWriteFast(PIN_STIM_ON, LOW);

  // --- Serial ---
  Serial.begin(115200);

  // --- Compute DDS phase scale factor ---
  phaseScale = 4294967296.0f / actualPWMFreq;

  // --- Generate linear gamma LUTs (default) ---
  generateLinearGammaLUT(gammaLUT1, PWM_MAX);
  generateLinearGammaLUT(gammaLUT2, PWM_MAX);
  currentGammaLUT_A = gammaLUT1;
  currentGammaLUT_B = gammaLUT1;

  // --- Generate raised cosine envelope LUT (0 to FIXED_SCALE) ---
  for (int i = 0; i < FADE_LUT_SIZE; i++) {
    float angle = PI * (float)i / (float)FADE_LUT_MAX;
    float val = (1.0f - cosf(angle)) / 2.0f;
    raisedCosineLUT[i] = (uint16_t)(val * (float)FIXED_SCALE);
    if (raisedCosineLUT[i] > FIXED_SCALE) raisedCosineLUT[i] = FIXED_SCALE;
  }

  // --- Set brightness range ---
  maxBrightness = PWM_MAX;
  midBrightness = PWM_MAX / 2;

  // --- Generate sine wave LUT ---
  generateSineWaveTable(maxBrightness);

  // --- Set default 50% duty cycle ---
  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

// ========================== MAIN LOOP ==========================

void loop() {
  GetSerialInput();
  if (newData) {
    newData = false;
    ActionSerial();
  }
}

// ========================== SERIAL INPUT ==========================

void GetSerialInput() {
  static byte ndx = 0;
  char endMarker = '\r';
  char rc;
  if (Serial.available() > 0) {
    rc = Serial.read();
    if (rc != endMarker) {
      receivedChars[ndx] = rc;
      ndx++;
      if (ndx >= numChars) {
        ndx = numChars - 1;
      }
    } else {
      receivedChars[ndx] = '\0';
      ndx = 0;
      newData = true;
    }
  }
}

// ========================== SERIAL COMMAND DISPATCH ==========================

void ActionSerial() {
  Serial.print(F("rc: "));
  Serial.print(receivedChars);
  Serial.print(F("\n"));

  char *token;
  uint8_t idx = 0;
#define MAX_VALS 50
  char *serialVals[MAX_VALS];
  token = strtok(receivedChars, ",");

  while (token != NULL) {
    if (idx < MAX_VALS)
      serialVals[idx++] = token;
    token = strtok(NULL, ",");
  }

  // Guard: if no tokens were parsed (empty/whitespace-only command), do nothing
  if (idx == 0) {
    memset(receivedChars, '\0', sizeof(receivedChars));
    return;
  }

  char *FirstChar = serialVals[0];

  // s, durationMs, freqA, freqB, phaseA, phaseB, contrastA, contrastB  [8 tokens]
  if (strcmp(FirstChar, "s") == 0)
  {
    if (idx < 8) { Serial.print(F("ERR: s needs 7 params\n")); return; }
    long stimulusDuration = atof(serialVals[1]);
    float freqA = atof(serialVals[2]);
    float freqB = atof(serialVals[3]);
    float phaseA = atof(serialVals[4]);
    float phaseB = atof(serialVals[5]);
    float contrastA = atof(serialVals[6]);
    float contrastB = atof(serialVals[7]);

    outputSinewave(freqA, freqB, stimulusDuration, phaseA, phaseB, contrastA, contrastB);
  }
  // se, durationMs, freqA, freqB, envFreq, maxContrastA, maxContrastB  [7 tokens]
  else if (strcmp(FirstChar, "se") == 0)
  {
    if (idx < 7) { Serial.print(F("ERR: se needs 6 params\n")); return; }
    long stimulusDuration = atof(serialVals[1]);
    float freqA = atof(serialVals[2]);
    float freqB = atof(serialVals[3]);
    float envFrequency = atof(serialVals[4]);
    float maxContrastA = atof(serialVals[5]);
    float maxContrastB = atof(serialVals[6]);

    SineContrastConv(stimulusDuration, freqA, freqB, envFrequency, maxContrastA, maxContrastB);
  }
  // fs, fmin, fmax, sweep, phA, phB, cA, cB  [8 tokens]
  else if (strcmp(FirstChar, "fs") == 0)
  {
    if (idx < 8) { Serial.print(F("ERR: fs needs 7 params\n")); return; }
    float fmin = atof(serialVals[1]);
    float fmax = atof(serialVals[2]);
    float sweepFactorPerSec = atof(serialVals[3]);
    float phaseA = atof(serialVals[4]);
    float phaseB = atof(serialVals[5]);
    float contrastA = atof(serialVals[6]);
    float contrastB = atof(serialVals[7]);

    FrequencySweep(fmin, fmax, sweepFactorPerSec,
                   phaseA, phaseB, contrastA, contrastB);
  }
  // sfs, f0, f1, step, cyc, phA, phB, cA, cB  [9 tokens]
  else if (strcmp(FirstChar, "sfs") == 0)
  {
    if (idx < 9) { Serial.print(F("ERR: sfs needs 8 params\n")); return; }
    float startFreq = atof(serialVals[1]);
    float endFreq = atof(serialVals[2]);
    float stepFreq = atof(serialVals[3]);
    int cyclesPerFreq = atoi(serialVals[4]);
    float phaseA = atof(serialVals[5]);
    float phaseB = atof(serialVals[6]);
    float contrastA = atof(serialVals[7]);
    float contrastB = atof(serialVals[8]);

    SteppedFrequencySweep(startFreq, endFreq, stepFreq, cyclesPerFreq, phaseA, phaseB, contrastA, contrastB);
  }
  // sd, dutyA, dutyB  [3 tokens]
  else if (strcmp(FirstChar, "sd") == 0)
  {
    if (idx < 3) { Serial.print(F("ERR: sd needs 2 params\n")); return; }
    float dutyCycle_A = atof(serialVals[1]);
    float dutyCycle_B = atof(serialVals[2]);
    setDutyCycle(dutyCycle_A, dutyCycle_B);
  }
  // sdt, dutyA, dutyB, durationMs  [4 tokens]
  else if (strcmp(FirstChar, "sdt") == 0)
  {
    if (idx < 4) { Serial.print(F("ERR: sdt needs 3 params\n")); return; }
    float dutyCycle_A = atof(serialVals[1]);
    float dutyCycle_B = atof(serialVals[2]);
    long stimulusDuration = atof(serialVals[3]);

    setDutyCycleTime(dutyCycle_A, dutyCycle_B, stimulusDuration);
  }
  // gc, step, wait, nReps  [4 tokens]
  else if (strcmp(FirstChar, "gc") == 0)
  {
    if (idx < 4) { Serial.print(F("ERR: gc needs 3 params\n")); return; }
    float stepSize = atof(serialVals[1]);
    long waitTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    cycleDutyCycles(stepSize, waitTime, nReps);
  }
  else if (strcmp(FirstChar, "useChB") == 0)
  {
    useChB = atoi(serialVals[1]);
    if (!useChB) {
      Serial.print(F("ChB OFF\n"));
    } else {
      Serial.print(F("ChB ON\n"));
    }
  }
  else if (strcmp(FirstChar, "useChA") == 0)
  {
    useChA = atoi(serialVals[1]);
    if (!useChA) {
      Serial.print(F("ChA OFF\n"));
    } else {
      Serial.print(F("ChA ON\n"));
    }
  }
  else if (strcmp(FirstChar, "ana") == 0) {
    readAnalogVals();
  }
  else if (strcmp(FirstChar, "agc") == 0)
  {
    uint8_t lutIndex = atoi(serialVals[1]);
    if (lutIndex == 1) {
      currentGammaLUT_A = gammaLUT1;
      currentGammaLUT_B = gammaLUT1;
      Serial.print(F("LUT 1 SELECTED\n"));
    } else if (lutIndex == 2) {
      currentGammaLUT_A = gammaLUT2;
      currentGammaLUT_B = gammaLUT2;
      Serial.print(F("LUT 2 SELECTED\n"));
    }
  }
  // bright,A/B,percent  — immediately rescales current PWM output
  else if (strcmp(FirstChar, "bright") == 0)
  {
    if (idx < 3) { Serial.print(F("ERR: bright needs 2 params\n")); return; }
    char *channel = serialVals[1];
    float percent = atof(serialVals[2]);
    if (percent < 0.0f) percent = 0.0f;
    if (percent > 100.0f) percent = 100.0f;
    uint16_t scale = (uint16_t)(percent / 100.0f * FIXED_SCALE);

    if (strcmp(channel, "A") == 0 || strcmp(channel, "a") == 0) {
      brightnessScaleA = scale;
      setChA(lastPwmValA);  // Immediately re-apply current output at new brightness
      Serial.print(F("ChA brightness: "));
    } else {
      brightnessScaleB = scale;
      setChB(lastPwmValB);  // Immediately re-apply current output at new brightness
      Serial.print(F("ChB brightness: "));
    }
    Serial.print(percent);
    Serial.print(F("%\n"));
  }
  // --- NEW: Upload empirical gamma LUT via ASCII ---
  // loadlut,A,v0,v1,v2,...,v4095
  else if (strcmp(FirstChar, "loadlut") == 0)
  {
    // Loads empirical gamma LUT data into gammaLUT1 (channel A) or gammaLUT2 (channel B).
    // Always writes into the named backing array, not into wherever the pointer currently
    // points — so the active LUT is only updated once loading is complete.
    //
    // Usage:
    //   "loadlut,A"              → stream mode: send comma-separated values then "END"
    //   "loadlut,B"              → same for channel B
    //   "loadlut,A,0,1,3,7,..."  → inline mode (max ~40 values due to 50-char buffer;
    //                              use stream mode for the full 4096-entry LUT)
    char *channel = serialVals[1];
    bool isChA = (strcmp(channel, "A") == 0 || strcmp(channel, "a") == 0);
    uint16_t *targetLUT = isChA ? gammaLUT1 : gammaLUT2;

    if (idx > 2) {
      // Inline mode — limited by numChars buffer, useful for small test LUTs
      int count = idx - 2;
      if (count > PWM_MAX + 1) count = PWM_MAX + 1;
      for (int i = 0; i < count; i++) {
        targetLUT[i] = (uint16_t)atoi(serialVals[i + 2]);
      }
      // Switch active pointer to the newly written slot
      if (isChA) { currentGammaLUT_A = targetLUT; }
      else        { currentGammaLUT_B = targetLUT; }
      Serial.print(F("LUT loaded: "));
      Serial.print(count);
      Serial.print(F(" values\n"));
    } else {
      // Stream mode: read 4096 comma-separated values then "END"
      Serial.print(F("SEND_LUT_DATA\n"));
      loadLUTFromSerial(targetLUT);
      // Switch active pointer only after successful load
      if (isChA) { currentGammaLUT_A = targetLUT; }
      else        { currentGammaLUT_B = targetLUT; }
    }
  }
  else
  {
    Serial.print(FirstChar);
    Serial.print(F(" is an invalid stimulus code - make sure you are using carriage return line ending\n"));
  }
  memset(receivedChars, '\0', sizeof(receivedChars));
}

// ========================== LUT UPLOAD (STREAM MODE) ==========================
// Reads 4096 comma-separated values from serial, split across multiple lines.
// Send values separated by commas, terminated with a line containing "END".

void loadLUTFromSerial(uint16_t *targetLUT) {
  int loaded = 0;
  char buf[64];
  int bufIdx = 0;
  unsigned long lastRxTime = millis();  // Timeout tracker

  while (loaded < PWM_MAX + 1) {

    // Escape hatch: if no data arrives within 5 seconds, abort
    if (millis() - lastRxTime > 5000) {
      Serial.print(F("LUT_LOAD_ERROR: TIMEOUT after "));
      Serial.print(loaded);
      Serial.print(F(" values\n"));
      return;  // Leaves targetLUT partially written; caller will not update the active pointer
    }

    if (Serial.available() > 0) {
      lastRxTime = millis();  // Reset timeout on every received byte
      char c = Serial.read();

      if (c == ',' || c == '\r' || c == '\n') {
        if (bufIdx > 0) {
          buf[bufIdx] = '\0';

          // Check for END marker
          if (strcmp(buf, "END") == 0 || strcmp(buf, "end") == 0) {
            break;
          }

          targetLUT[loaded] = (uint16_t)atoi(buf);
          loaded++;
          bufIdx = 0;
        }
      } else if (bufIdx < 63) {
        buf[bufIdx++] = c;
      }
    }
  }

  Serial.print(F("LUT_LOADED: "));
  Serial.print(loaded);
  Serial.print(F(" values\n"));
}

// ========================== GAMMA LUT GENERATION ==========================

void generateLinearGammaLUT(uint16_t *lut, uint16_t maxOutput) {
  for (int i = 0; i <= PWM_MAX; i++) {
    lut[i] = (uint16_t)((uint32_t)i * maxOutput / PWM_MAX);
  }
}

void generateGammaLUT(uint16_t *lut, float gamma, uint16_t maxOutput) {
  for (int i = 0; i <= PWM_MAX; i++) {
    float normalized = (float)i / (float)PWM_MAX;
    float corrected = powf(normalized, gamma);
    lut[i] = (uint16_t)(corrected * maxOutput + 0.5f);
  }
}

// ========================== SINE WAVETABLE GENERATION ==========================

void generateSineWaveTable(uint16_t top) {
  for (int i = 0; i < TABLE_SIZE; i++) {
    float angle = (2.0f * PI * i) / TABLE_SIZE;
    // +0.5f for rounding to nearest integer instead of truncating
    sineWaveTable[i] = (uint16_t)(((sinf(angle) + 1.0f) * (top / 2.0f)) + 0.5f);
  }
}

// ========================== SINEWAVE FLICKER (DDS) ==========================

void outputSinewave(float freqA, float freqB, long duration, float phaseA, float phaseB, float contrastA, float contrastB) {

  // Phase: 0.0 = 0°, 0.25 = 90°, 0.5 = 180°, etc.
  phaseAccumulatorA = (uint32_t)(phaseA * 4294967296.0f);
  phaseAccumulatorB = (uint32_t)(phaseB * 4294967296.0f);

  contrastMultIntA = (uint16_t)(contrastA * (float)FIXED_SCALE);
  contrastMultIntB = (uint16_t)(contrastB * (float)FIXED_SCALE);

  phaseIncrementA = calcPhaseInc(freqA);
  phaseIncrementB = calcPhaseInc(freqB);

  // --- Temporal window & master clock setup ---
  unsigned long totalInterrupts = (unsigned long)((duration / 1000.0f) * actualPWMFreq);

  float fadeTimeMs = FADE_DURATION_MS;

  if (fadeTimeMs > 0.0f) {
    fadeInterrupts = (unsigned long)(actualPWMFreq * (fadeTimeMs / 1000.0f));
    if (fadeInterrupts > totalInterrupts / 2) {
      fadeInterrupts = totalInterrupts / 2;
    }
    fadeOutStartInterrupt = totalInterrupts - fadeInterrupts;

    // Mathematically forces division to round UP
    envStep = (((uint32_t)FADE_LUT_MAX << 16) + fadeInterrupts - 1) / fadeInterrupts;
  } else {
    fadeInterrupts = 0;
    fadeOutStartInterrupt = 4294967295UL;
    envStep = 0;
  }

  currentTick = 0;
  completedCycles = 0;

  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }

  digitalWriteFast(PIN_STIM_ON, HIGH);
  digitalWriteFast(PIN_CYCLE_INDICATOR, HIGH);

  startDDS(sinewaveInterrupt);

  // Wait for the master clock to reach total duration
  unsigned long safeTick = 0;
  while (true) {
    noInterrupts();
    safeTick = currentTick;
    interrupts();

    if (safeTick >= totalInterrupts) break;

    delayMicroseconds(1);
  }

  stopDDS();

  digitalWriteFast(PIN_CYCLE_INDICATOR, LOW);
  digitalWriteFast(PIN_STIM_ON, LOW);
  Serial.print(F("-1\n"));

  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

void sinewaveInterrupt() {
  uint16_t currentEnvelope = FIXED_SCALE;

  if (currentTick < fadeInterrupts) {
    uint32_t phase = currentTick * envStep;
    uint16_t envIndex = phase >> 16;
    if (envIndex > FADE_LUT_MAX) envIndex = FADE_LUT_MAX;
    currentEnvelope = raisedCosineLUT[envIndex];

  } else if (currentTick >= fadeOutStartInterrupt) {
    uint32_t ticksIntoFadeOut = currentTick - fadeOutStartInterrupt;
    uint32_t phase = ticksIntoFadeOut * envStep;
    uint16_t shiftPhase = phase >> 16;
    uint16_t envIndex = (shiftPhase >= FADE_LUT_MAX) ? 0 : (FADE_LUT_MAX - shiftPhase);
    currentEnvelope = raisedCosineLUT[envIndex];
  }
  currentTick++;

  uint16_t effectiveContrastA = ((uint32_t)contrastMultIntA * currentEnvelope) >> FIXED_SHIFT;
  uint16_t effectiveContrastB = ((uint32_t)contrastMultIntB * currentEnvelope) >> FIXED_SHIFT;

  uint32_t oldPhaseA = phaseAccumulatorA;

  phaseAccumulatorA += phaseIncrementA;
  phaseAccumulatorB += phaseIncrementB;

  // Cycle indicator toggle on Channel A phase wraparound
  if (phaseAccumulatorA < oldPhaseA) {
    completedCycles++;
    digitalToggleFast(PIN_CYCLE_INDICATOR);
  }

  uint16_t indexA = phaseAccumulatorA >> 20;  // Top 12 bits → 4096 entries
  uint16_t indexB = phaseAccumulatorB >> 20;

  int32_t tempA = (int32_t)sineWaveTable[indexA] - midBrightness;
  int32_t pwmValA_calc = midBrightness + ((tempA * effectiveContrastA) >> FIXED_SHIFT);
  if (pwmValA_calc < 0) pwmValA_calc = 0;
  if (pwmValA_calc > PWM_MAX) pwmValA_calc = PWM_MAX;

  int32_t tempB = (int32_t)sineWaveTable[indexB] - midBrightness;
  int32_t pwmValB_calc = midBrightness + ((tempB * effectiveContrastB) >> FIXED_SHIFT);
  if (pwmValB_calc < 0) pwmValB_calc = 0;
  if (pwmValB_calc > PWM_MAX) pwmValB_calc = PWM_MAX;

  if (useChA) { setChA((uint16_t)pwmValA_calc); }
  if (useChB) { setChB((uint16_t)pwmValB_calc); }
}

// ========================== SINE WITH CONTRAST ENVELOPE (DDS) ==========================

void SineContrastConv(float duration, float freqA, float freqB, float envelopeFreq, float maxContrastA, float maxContrastB) {

  contrastMultIntA = (uint16_t)(maxContrastA * (float)FIXED_SCALE);
  contrastMultIntB = (uint16_t)(maxContrastB * (float)FIXED_SCALE);

  phaseAccumulatorA = 0;
  phaseAccumulatorB = 0;
  phaseIncrementA = calcPhaseInc(freqA);
  phaseIncrementB = calcPhaseInc(freqB);

  envAccumulator = (192UL << 24);  // Start at 270° (absolute trough)
  envIncrement = calcPhaseInc(envelopeFreq);

  fadeInterrupts = 0;
  fadeOutStartInterrupt = 4294967295UL;

  unsigned long totalInterrupts = (unsigned long)((duration / 1000.0f) * actualPWMFreq);
  currentTick = 0;
  completedCycles = 0;

  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }

  digitalWriteFast(PIN_STIM_ON, HIGH);
  digitalWriteFast(PIN_CYCLE_INDICATOR, HIGH);  // Raise indicator to match all other stim functions

  startDDS(sinewaveEnvelopeInterrupt);

  unsigned long safeTick = 0;
  while (true) {
    noInterrupts();
    safeTick = currentTick;
    interrupts();

    if (safeTick >= totalInterrupts) break;
    delayMicroseconds(1);
  }

  stopDDS();

  digitalWriteFast(PIN_CYCLE_INDICATOR, LOW);
  digitalWriteFast(PIN_STIM_ON, LOW);
  Serial.print(F("-1\n"));
  Serial.flush();

  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

void sinewaveEnvelopeInterrupt() {
  currentTick++;

  uint32_t oldPhaseA = phaseAccumulatorA;

  phaseAccumulatorA += phaseIncrementA;
  phaseAccumulatorB += phaseIncrementB;

  if (phaseAccumulatorA < oldPhaseA) {
    completedCycles++;
  }

  uint32_t oldEnvAccumulator = envAccumulator;
  envAccumulator += envIncrement;

  if (envAccumulator < oldEnvAccumulator) {
    digitalToggleFast(PIN_CYCLE_INDICATOR);
  }

  uint16_t indexA = phaseAccumulatorA >> 20;
  uint16_t indexB = phaseAccumulatorB >> 20;
  uint16_t indexEnv = envAccumulator >> 20;

  // Use sine LUT for envelope — gives raised-cosine-like contrast modulation
  uint16_t contrastMultInt = sineWaveTable[indexEnv];

  uint16_t currentContrastIntA = ((uint32_t)contrastMultInt * contrastMultIntA) >> FIXED_SHIFT;
  uint16_t currentContrastIntB = ((uint32_t)contrastMultInt * contrastMultIntB) >> FIXED_SHIFT;

  int32_t tempA = (int32_t)sineWaveTable[indexA] - midBrightness;
  int32_t pwmValA_calc = midBrightness + ((tempA * currentContrastIntA) >> FIXED_SHIFT);
  if (pwmValA_calc < 0) pwmValA_calc = 0;
  if (pwmValA_calc > PWM_MAX) pwmValA_calc = PWM_MAX;

  int32_t tempB = (int32_t)sineWaveTable[indexB] - midBrightness;
  int32_t pwmValB_calc = midBrightness + ((tempB * currentContrastIntB) >> FIXED_SHIFT);
  if (pwmValB_calc < 0) pwmValB_calc = 0;
  if (pwmValB_calc > PWM_MAX) pwmValB_calc = PWM_MAX;

  if (useChA) { setChA((uint16_t)pwmValA_calc); }
  if (useChB) { setChB((uint16_t)pwmValB_calc); }
}

// ========================== FREQUENCY SWEEP (STEPPED) ==========================

void SteppedFrequencySweep(float startFreq, float endFreq, float stepFreq, int cyclesPerFreq, float phaseA, float phaseB, float contrastA, float contrastB) {

  phaseAccumulatorA = (uint32_t)(phaseA * 4294967296.0f);
  phaseAccumulatorB = (uint32_t)(phaseB * 4294967296.0f);
  contrastMultIntA = (uint16_t)(contrastA * (float)FIXED_SCALE);
  contrastMultIntB = (uint16_t)(contrastB * (float)FIXED_SCALE);

  fadeInterrupts = 0;
  fadeOutStartInterrupt = 4294967295UL;

  digitalWriteFast(PIN_STIM_ON, HIGH);
  digitalWriteFast(PIN_CYCLE_INDICATOR, HIGH);

  float currentFreq = startFreq;
  bool sweepingUp = startFreq <= endFreq;
  stepFreq = fabsf(stepFreq);

  phaseIncrementA = calcPhaseInc(currentFreq);
  phaseIncrementB = phaseIncrementA;

  startDDS(sinewaveInterrupt);

  while ((sweepingUp && currentFreq <= endFreq) || (!sweepingUp && currentFreq >= endFreq)) {

    uint32_t newInc = calcPhaseInc(currentFreq);
    noInterrupts();
    phaseIncrementA = newInc;
    phaseIncrementB = newInc;
    interrupts();

    uint32_t startCycles = 0;
    noInterrupts(); startCycles = completedCycles; interrupts();

    uint32_t safeCycles = 0;
    while (true) {
      noInterrupts();
      safeCycles = completedCycles;
      interrupts();

      if ((safeCycles - startCycles) >= (uint32_t)cyclesPerFreq) break;
      delayMicroseconds(1);
    }

    if (sweepingUp) {
      currentFreq += stepFreq;
    } else {
      currentFreq -= stepFreq;
    }
  }

  stopDDS();

  digitalWriteFast(PIN_CYCLE_INDICATOR, LOW);
  digitalWriteFast(PIN_STIM_ON, LOW);
  Serial.print(F("-1\n"));

  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

// ========================== FREQUENCY SWEEP (CONTINUOUS) ==========================

void FrequencySweep(float fmin, float fmax, float sweepFactorPerSec,
                    float phaseA, float phaseB, float contrastA, float contrastB) {

  if (fmin <= 0.0f) fmin = 0.001f;
  if (fmax <= 0.0f) fmax = 0.001f;

  if (fmin > fmax) {
    float temp = fmin;
    fmin = fmax;
    fmax = temp;
  }

  phaseAccumulatorA = (uint32_t)(phaseA * 4294967296.0f);
  phaseAccumulatorB = (uint32_t)(phaseB * 4294967296.0f);
  contrastMultIntA = (uint16_t)(contrastA * (float)FIXED_SCALE);
  contrastMultIntB = (uint16_t)(contrastB * (float)FIXED_SCALE);

  fadeInterrupts = 0;
  fadeOutStartInterrupt = 4294967295UL;

  unsigned long totalDurationUs = 0;
  bool isSweeping = false;

  if (sweepFactorPerSec > 0.000001f && fmax > fmin && (fmax / fmin) > 1.000001f) {
    float timeForOneWaySweepSec_calc = logf(fmax / fmin) / sweepFactorPerSec;

    if (timeForOneWaySweepSec_calc > 0.0000001f) {
      isSweeping = true;
      totalDurationUs = 2UL * (unsigned long)(timeForOneWaySweepSec_calc * 1000000.0f);
    }
  }

  if (!isSweeping || totalDurationUs == 0) {
    totalDurationUs = TARGET_RECONFIG_INTERVAL_US;
    isSweeping = false;
  }

  float actual_calc_freq = fmin;

  phaseIncrementA = calcPhaseInc(actual_calc_freq);
  phaseIncrementB = phaseIncrementA;

  digitalWriteFast(PIN_STIM_ON, HIGH);
  digitalWriteFast(PIN_CYCLE_INDICATOR, HIGH);

  unsigned long totalInterrupts = (unsigned long)((totalDurationUs / 1000000.0f) * actualPWMFreq);
  currentTick = 0;
  completedCycles = 0;

  startDDS(sinewaveInterrupt);

  unsigned long safeTick = 0;

  while (true) {
    noInterrupts();
    safeTick = currentTick;
    interrupts();

    if (safeTick >= totalInterrupts) break;

    float sweepProgress = (float)safeTick / (float)totalInterrupts;

    if (isSweeping) {
      if (sweepProgress < 0.5f) {
        actual_calc_freq = fmin * expf(sweepFactorPerSec * (sweepProgress * totalDurationUs / 1000000.0f));
        if (actual_calc_freq > fmax) actual_calc_freq = fmax;
      } else {
        actual_calc_freq = fmax * expf(-sweepFactorPerSec * ((sweepProgress - 0.5f) * totalDurationUs / 1000000.0f));
        if (actual_calc_freq < fmin) actual_calc_freq = fmin;
      }

      uint32_t newInc = calcPhaseInc(actual_calc_freq);
      noInterrupts();
      phaseIncrementA = newInc;
      phaseIncrementB = newInc;
      interrupts();
    }

    delayMicroseconds(TARGET_RECONFIG_INTERVAL_US);
  }

  stopDDS();

  digitalWriteFast(PIN_CYCLE_INDICATOR, LOW);
  digitalWriteFast(PIN_STIM_ON, LOW);
  Serial.print(F("-1\n"));

  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

// ========================== STATIC DUTY CYCLE FUNCTIONS ==========================

void setDutyCycle(float dutyCyclePercentage_A, float dutyCyclePercentage_B) {
  if (dutyCyclePercentage_A < 0.0f) dutyCyclePercentage_A = 0.0f;
  if (dutyCyclePercentage_A > 100.0f) dutyCyclePercentage_A = 100.0f;
  uint16_t pwmValueA = (uint16_t)((dutyCyclePercentage_A / 100.0f) * maxBrightness);

  if (dutyCyclePercentage_B < 0.0f) dutyCyclePercentage_B = 0.0f;
  if (dutyCyclePercentage_B > 100.0f) dutyCyclePercentage_B = 100.0f;
  uint16_t pwmValueB = (uint16_t)((dutyCyclePercentage_B / 100.0f) * maxBrightness);

  digitalToggleFast(PIN_CYCLE_INDICATOR);
  if (useChA) { setChA(pwmValueA); }
  if (useChB) { setChB(pwmValueB); }
}

void setDutyCycleTime(float dutyCyclePercentage_A, float dutyCyclePercentage_B, long duration) {
  if (dutyCyclePercentage_A < 0.0f) dutyCyclePercentage_A = 0.0f;
  if (dutyCyclePercentage_A > 100.0f) dutyCyclePercentage_A = 100.0f;
  uint16_t pwmValueA = (uint16_t)((dutyCyclePercentage_A / 100.0f) * maxBrightness);

  if (dutyCyclePercentage_B < 0.0f) dutyCyclePercentage_B = 0.0f;
  if (dutyCyclePercentage_B > 100.0f) dutyCyclePercentage_B = 100.0f;
  uint16_t pwmValueB = (uint16_t)((dutyCyclePercentage_B / 100.0f) * maxBrightness);

  unsigned long startTime = millis();  // unsigned long to avoid signed overflow after ~25 days
  digitalWriteFast(PIN_STIM_ON, HIGH);
  digitalWriteFast(PIN_CYCLE_INDICATOR, HIGH);
  if (useChA) { setChA(pwmValueA); }
  if (useChB) { setChB(pwmValueB); }
  while (millis() - startTime < (unsigned long)duration) {
    delayMicroseconds(1);
  }
  digitalWriteFast(PIN_CYCLE_INDICATOR, LOW);
  digitalWriteFast(PIN_STIM_ON, LOW);
  Serial.print(F("-1\n"));
  Serial.flush();
  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

void cycleDutyCycles(float stepSize, float waitTime, int nReps) {
  for (int irep = 0; irep < nReps; irep++) {
    float dutyCycle = 0;
    while (dutyCycle <= 1.0f) {
      Serial.print(dutyCycle);
      Serial.print(F("\n"));
      uint16_t pwmValue = (uint16_t)(dutyCycle * maxBrightness);
      if (useChA) { setChA(pwmValue); }
      if (useChB) { setChB(pwmValue); }
      dutyCycle = dutyCycle + stepSize;
      delay(waitTime);
    }
  }
  Serial.print(F("-1\n"));
  if (useChA) { setChA(midBrightness); }
  if (useChB) { setChB(midBrightness); }
}

// ========================== ANALOG READ ==========================

void readAnalogVals() {
  bool keepReading = true;
  const unsigned long interval = 100;
  unsigned long previousMillis = millis();

  // Teensy 4.1: A0 = pin 14, A1 = pin 15 (same A0/A1 names work)
  int analogValue0 = analogRead(A0);
  int analogValue1 = analogRead(A1);

  while (keepReading) {
    unsigned long currentMillis = millis();
    if (currentMillis - previousMillis >= interval) {
      previousMillis = currentMillis;

      analogValue0 = analogRead(A0);
      analogValue1 = analogRead(A1);
      Serial.print(analogValue0);
      Serial.print(F(","));
      Serial.println(analogValue1);
    }

    if (Serial.available() > 0) {
      String input = Serial.readStringUntil('\n');
      input.trim();

      if (input.equalsIgnoreCase("done")) {
        keepReading = false;
        Serial.println(F("Stopped reading analog values."));
      }
      // Dynamic LUT switching based on analog reading (preserved from Leonardo)
      if (analogValue0 < 420) {
        currentGammaLUT_A = gammaLUT2;
        currentGammaLUT_B = gammaLUT2;
        Serial.print(F("LUT 2 SELECTED\n"));
      } else {
        currentGammaLUT_A = gammaLUT1;
        currentGammaLUT_B = gammaLUT1;
        Serial.print(F("LUT 1 SELECTED\n"));
      }
    }
  }
}

// ========================== INTERVAL TIMER CONTROL ==========================

void startDDS(void (*isrFunc)()) {
  float intervalUs = 1000000.0f / actualPWMFreq;
  // begin() must be called BEFORE priority() — begin() allocates the hardware PIT channel
  // and sets up the NVIC vector. priority() then updates that vector's priority level.
  // Priority 16 preempts Serial, millis, and most other system interrupts (default 128).
  // Lower number = higher priority. 0 is highest, 255 is lowest.
  ddsTimer.begin(isrFunc, intervalUs);
  ddsTimer.priority(16);
}

void stopDDS() {
  ddsTimer.end();
}
