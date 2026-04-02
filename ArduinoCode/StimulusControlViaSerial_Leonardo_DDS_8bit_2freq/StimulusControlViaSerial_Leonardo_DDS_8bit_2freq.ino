#include <math.h>  

#define CLOCK_FREQ 16000000                // Arduino Leonardo clock frequency (16 MHz)
#define TABLE_SIZE 256                     // Number of samples in the wavetable
#define TARGET_RECONFIG_INTERVAL_US 1000L  // FOR FREQUENCY SWEEP 1000 microseconds = 1 millisecond

// --- FADE & ENVELOPE MACROS ---
#define FADE_DURATION_MS 100.0             // Easily change the fade-in/out time in milliseconds here!
#define FADE_LUT_SIZE 64                   // Change this to 64, 128, 256, 512, etc.
#define FADE_LUT_MAX (FADE_LUT_SIZE - 1)   // Used for safe zero-indexed math

// 8-Bit Gamma-correction LUTS (Replace with your calibrated 256-item arrays!)
const uint16_t PROGMEM ChA1LUT[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
  2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
  6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 11, 11,
  11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19,
  19, 20, 21, 21, 22, 22, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29,
  30, 31, 31, 32, 33, 34, 34, 35, 36, 37, 38, 38, 39, 40, 41, 42,
  43, 44, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
  58, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 73, 74, 75,
  76, 77, 79, 80, 81, 82, 84, 85, 86, 88, 89, 90, 92, 93, 94, 96,
  97, 98, 100, 101, 103, 104, 106, 107, 109, 110, 112, 113, 115, 116, 118, 119,
  121, 123, 124, 126, 128, 129, 131, 133, 134, 136, 138, 139, 141, 143, 145, 147,
  148, 150, 152, 154, 156, 158, 160, 161, 163, 165, 167, 169, 171, 173, 175, 177,
  179, 181, 183, 185, 188, 190, 192, 194, 196, 198, 201, 203, 205, 207, 209, 212,
  214, 216, 218, 221, 223, 225, 228, 230, 233, 235, 237, 240, 242, 245, 247, 250,
  252, 255
};

const uint16_t PROGMEM ChA2LUT[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
  2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
  6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 11, 11,
  11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19,
  19, 20, 21, 21, 22, 22, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29,
  30, 31, 31, 32, 33, 34, 34, 35, 36, 37, 38, 38, 39, 40, 41, 42,
  43, 44, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
  58, 59, 60, 61, 62, 63, 65, 66, 67, 68, 69, 70, 71, 73, 74, 75,
  76, 77, 79, 80, 81, 82, 84, 85, 86, 88, 89, 90, 92, 93, 94, 96,
  97, 98, 100, 101, 103, 104, 106, 107, 109, 110, 112, 113, 115, 116, 118, 119,
  121, 123, 124, 126, 128, 129, 131, 133, 134, 136, 138, 139, 141, 143, 145, 147,
  148, 150, 152, 154, 156, 158, 160, 161, 163, 165, 167, 169, 171, 173, 175, 177,
  179, 181, 183, 185, 188, 190, 192, 194, 196, 198, 201, 203, 205, 207, 209, 212,
  214, 216, 218, 221, 223, 225, 228, 230, 233, 235, 237, 240, 242, 245, 247, 250,
  252, 255
};

const uint16_t *currentChALUT;
const uint16_t *currentChBLUT;

// PWM variables
long prescaler;
uint16_t TOP;  
long TopLumi;  
long MidLumi;

// 31.25 kHz target to kill saccadic CFF artifacts
long desiredPWMFrequency = 31250;  
float actualPWMFreq = 31250.0;     // Hardware-accurate frequency for DDS calculations

// channel selection
bool useChA = true;
bool useChB = true;

// serial
const byte numChars = 50;
char receivedChars[numChars];  
bool newData = false;

// Array to hold the wavetable
uint16_t sineWaveTable[TABLE_SIZE];

// --- DDS (Direct Digital Synthesis) Variables ---
volatile uint32_t phaseAccumulatorA = 0;
volatile uint32_t phaseAccumulatorB = 0;
volatile uint32_t phaseIncrementA = 0;
volatile uint32_t phaseIncrementB = 0;

volatile uint32_t envAccumulator = 0;
volatile uint32_t envIncrement = 0;

volatile uint16_t contrastMultIntA = 256;
volatile uint16_t contrastMultIntB = 256;

// Tracks full sine wave cycles safely up to 4.2 billion
volatile uint32_t completedCycles = 0; 

// --- Raised Cosine Temporal Window Variables ---
uint16_t raisedCosineLUT[FADE_LUT_SIZE]; 
volatile unsigned long currentTick = 0; // The ultimate Master Clock
volatile unsigned long fadeInterrupts = 0;
volatile unsigned long fadeOutStartInterrupt = 0;
volatile uint32_t envStep = 0; 

// Function to calculate exact DDS phase increment 
uint32_t calcPhaseInc(float freq) {
  uint64_t scaledFreq = (uint64_t)(freq * 10000.0); 
  return (uint32_t)((scaledFreq * 4294967296ULL) / ((uint64_t)(actualPWMFreq * 10000.0)));
}

void setup() {
  pinMode(9, OUTPUT);   // Pin 9 controlled by Timer1 (Channel A)
  pinMode(10, OUTPUT);  // Pin 10 controlled by Timer1 (Channel B)

  pinMode(4, OUTPUT);  // Pin 4 indicator pin
  pinMode(5, OUTPUT);  // Pin 5 stim ON or OFF pin

  PORTD &= ~(1 << PIND4);   
  PORTC &= ~(1 << PORTC6);  

  Serial.begin(115200);

  currentChALUT = ChA1LUT;
  currentChBLUT = ChA1LUT;

  // Generate the Modular Raised Cosine Envelope LUT (0 to 256 scale)
  for (int i = 0; i < FADE_LUT_SIZE; i++) {
    float angle = PI * (float)i / (float)FADE_LUT_MAX;         
    float val = (1.0 - cos(angle)) / 2.0;        
    raisedCosineLUT[i] = (uint16_t)(val * 256.0);
    if (raisedCosineLUT[i] > 256) raisedCosineLUT[i] = 256; 
  }

  // Calculate the prescaler and TOP value
  TOP = calculatePrescalerAndTOP(desiredPWMFrequency, prescaler);

  // Store the exact mathematical frequency the hardware actually achieved for DDS math
  actualPWMFreq = (float)CLOCK_FREQ / (2.0 * (float)prescaler * (float)TOP);

  // Apply the prescaler and TOP value to Timer1
  configureTimer1(prescaler, TOP);

  TopLumi = TOP;  
  MidLumi = TOP / 2;

  // Generate the sine wave LUT
  generateSineWaveTable(TopLumi);

  if (useChA) { setChA(TopLumi / 2); }  
  if (useChB) { setChB(TopLumi / 2); }  
}

void loop() {
  GetSerialInput();
  if (newData) {
    newData = false;
    ActionSerial();
  }
}

//////////////////////////////////////// HANDLE SERIAL INPUT //////////////////////////////////////
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

void ActionSerial() {  
  Serial.print(F("rc: "));
  Serial.print(receivedChars);
  Serial.print(F("\n"));
  char delimiters[] = ",";
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

  char *FirstChar = serialVals[0];

  // s, durationMs, freqA, freqB, phaseA, phaseB, contrastA, contrastB
  if (strcmp(FirstChar, "s") == 0)  
  {
    long stimulusDuration = atof(serialVals[1]);
    float freqA = atof(serialVals[2]);
    float freqB = atof(serialVals[3]);
    float phaseA = atof(serialVals[4]);
    float phaseB = atof(serialVals[5]);
    float contrastA = atof(serialVals[6]);
    float contrastB = atof(serialVals[7]);

    outputSinewave(freqA, freqB, stimulusDuration, phaseA, phaseB, contrastA, contrastB);

  } 
  // se, durationMs, freqA, freqB, envFreq, maxContrastA, maxContrastB
  else if (strcmp(FirstChar, "se") == 0)  
  {
    long stimulusDuration = atof(serialVals[1]);
    float freqA = atof(serialVals[2]);
    float freqB = atof(serialVals[3]);
    float envFrequency = atof(serialVals[4]);
    float maxContrastA = atof(serialVals[5]); 
    float maxContrastB = atof(serialVals[6]); 

    SineContrastConv(stimulusDuration, freqA, freqB, envFrequency, maxContrastA, maxContrastB);

  } else if (strcmp(FirstChar, "fs") == 0)  
  {
    float fmin = atof(serialVals[1]);
    float fmax = atof(serialVals[2]);
    float sweepFactorPerSec = atof(serialVals[3]);
    float phaseA = atof(serialVals[4]);
    float phaseB = atof(serialVals[5]);
    float contrastA = atof(serialVals[6]);
    float contrastB = atof(serialVals[7]);

    FrequencySweep(fmin, fmax, sweepFactorPerSec,
                   phaseA, phaseB, contrastA, contrastB);

  } else if (strcmp(FirstChar, "sfs") == 0)  
  {
    float startFreq = atof(serialVals[1]);
    float endFreq = atof(serialVals[2]);
    float stepFreq = atof(serialVals[3]);
    int cyclesPerFreq = atoi(serialVals[4]);
    float phaseA = atof(serialVals[5]);
    float phaseB = atof(serialVals[6]);
    float contrastA = atof(serialVals[7]);
    float contrastB = atof(serialVals[8]);

    SteppedFrequencySweep(startFreq, endFreq, stepFreq, cyclesPerFreq, phaseA, phaseB, contrastA, contrastB);

  } else if (strcmp(FirstChar, "sd") == 0)  
  {
    float dutyCycle_A = atof(serialVals[1]);
    float dutyCycle_B = atof(serialVals[2]);
    setDutyCycle(dutyCycle_A, dutyCycle_B, TopLumi);

  } else if (strcmp(FirstChar, "sdt") == 0)  
  {
    float dutyCycle_A = atof(serialVals[1]);
    float dutyCycle_B = atof(serialVals[2]);
    long stimulusDuration = atof(serialVals[3]);

    setDutyCycleTime(dutyCycle_A, dutyCycle_B, stimulusDuration, TopLumi);

  } else if (strcmp(FirstChar, "gc") == 0)  
  {
    float stepSize = atof(serialVals[1]);
    long waitTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    cycleDutyCycles(stepSize, waitTime, nReps, TopLumi);

  } else if (strcmp(FirstChar, "useChB") == 0)  
  {
    useChB = atoi(serialVals[1]);
    if (!useChB) {
      Serial.print(F("ChB OFF\n"));
    } else {
      Serial.print(F("ChB ON\n"));
    };

  } else if (strcmp(FirstChar, "useChA") == 0)  
  {
    useChA = atoi(serialVals[1]);
    if (!useChA) {
      Serial.print(F("ChA OFF\n"));
    } else {
      Serial.print(F("ChA ON\n"));
    };

  } else if (strcmp(FirstChar, "ana") == 0) {
    readAnalogVals();

  } else if (strcmp(FirstChar, "agc") == 0)  
  {
    uint8_t lutIndex = atoi(serialVals[1]);
    if (lutIndex == 1) {
      currentChALUT = ChA1LUT;
      currentChBLUT = ChA1LUT;
      Serial.print(F("LUT 1 SELECTED\n"));

    } else if (lutIndex == 2) {
      currentChALUT = ChA2LUT;
      currentChBLUT = ChA2LUT;
      Serial.print(F("LUT 2 SELECTED\n"));
    }
  } else  
  {
    Serial.print(FirstChar);
    Serial.print(F(" is an invalid stimulus code - make sure you are using carriage return line ending\n"));
  }
  memset(receivedChars, '\0', sizeof(receivedChars));
}

///////////////////////////////////// SINEWAVE FLICKER (DDS) //////////////////////////////////////
void generateSineWaveTable(long TOP) {
  for (int i = 0; i < TABLE_SIZE; i++) {
    float angle = (2.0 * PI * i) / TABLE_SIZE;                        
    sineWaveTable[i] = (uint16_t)((sin(angle) + 1.0) * (TOP / 2.0));  
  }
}

void outputSinewave(float freqA, float freqB, long duration, float phaseA, float phaseB, float contrastA, float contrastB) {

  phaseAccumulatorA = (uint32_t)(phaseA * 4294967296.0);
  phaseAccumulatorB = (uint32_t)(phaseB * 4294967296.0);

  contrastMultIntA = (uint16_t)(contrastA * 256.0);
  contrastMultIntB = (uint16_t)(contrastB * 256.0);

  phaseIncrementA = calcPhaseInc(freqA);
  phaseIncrementB = calcPhaseInc(freqB);

  // --- TEMPORAL WINDOW & MASTER CLOCK SETUP ---
  unsigned long totalInterrupts = (unsigned long)((duration / 1000.0) * actualPWMFreq); 
  
  float fadeTimeMs = FADE_DURATION_MS; 
  
  if (fadeTimeMs > 0.0) {
    fadeInterrupts = (unsigned long)(actualPWMFreq * (fadeTimeMs / 1000.0));
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
  
  currentTick = 0; // Master Clock Reset
  completedCycles = 0;
  // -------------------------------------------------

  if (useChA) { setChA(MidLumi); }
  if (useChB) { setChB(MidLumi); }

  PORTC |= (1 << PORTC6);                     
  PORTD |= (1 << PIND4);                      

  TIMSK0 &= ~_BV(TOIE0); 
  
  setTimer1Callback(sinewaveInterrupt);
  startTimer1Interrupt(); 

  // ATOMIC READ: Waiting for the Master Clock to reach total duration
  unsigned long safeTick = 0;
  while (true) {
    noInterrupts();
    safeTick = currentTick;
    interrupts();
    
    if (safeTick >= totalInterrupts) break;
    
    delayMicroseconds(1);  
  }

  stopTimer1Interrupt();    
  TIMSK0 |= _BV(TOIE0);     

  PORTD &= ~(1 << PIND4);   
  PORTC &= ~(1 << PORTC6);  
  Serial.print(F("-1\n"));

  if (useChA) { setChA(TopLumi / 2); }  
  if (useChB) { setChB(TopLumi / 2); }  
}

void sinewaveInterrupt() {
  uint16_t currentEnvelope = 256; 

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

  uint16_t effectiveContrastA = ((uint32_t)contrastMultIntA * currentEnvelope) >> 8;
  uint16_t effectiveContrastB = ((uint32_t)contrastMultIntB * currentEnvelope) >> 8;

  uint32_t oldPhaseA = phaseAccumulatorA;
  
  phaseAccumulatorA += phaseIncrementA;
  phaseAccumulatorB += phaseIncrementB;

  // Debug Pin toggle on Channel A cycles & cycle counting for sweeps
  if (phaseAccumulatorA < oldPhaseA) {
    completedCycles++;
    PORTD ^= (1 << PIND4); 
  }

  uint8_t indexA = phaseAccumulatorA >> 24;
  uint8_t indexB = phaseAccumulatorB >> 24;

  long tempA = (long)sineWaveTable[indexA] - MidLumi;
  long ocrValA_calc = MidLumi + ((tempA * effectiveContrastA) >> 8);
  if (ocrValA_calc < 0) ocrValA_calc = 0;
  if (ocrValA_calc > 255) ocrValA_calc = 255;
  uint16_t ocrValA = (uint16_t)ocrValA_calc; 
  
  long tempB = (long)sineWaveTable[indexB] - MidLumi;
  long ocrValB_calc = MidLumi + ((tempB * effectiveContrastB) >> 8);
  if (ocrValB_calc < 0) ocrValB_calc = 0;
  if (ocrValB_calc > 255) ocrValB_calc = 255;
  uint16_t ocrValB = (uint16_t)ocrValB_calc; 

  if (useChA) { setChA(ocrValA); }  
  if (useChB) { setChB(ocrValB); }  
}

/////////////////////////////////// SINE WAVE FLICKER WITH CONTRAST ENVELOPE (DDS) //////////////////////////
void SineContrastConv(float duration, float freqA, float freqB, float envelopeFreq, float maxContrastA, float maxContrastB) {

  contrastMultIntA = (uint16_t)(maxContrastA * 256.0);
  contrastMultIntB = (uint16_t)(maxContrastB * 256.0);

  phaseAccumulatorA = 0;
  phaseAccumulatorB = 0;
  phaseIncrementA = calcPhaseInc(freqA);
  phaseIncrementB = calcPhaseInc(freqB);

  envAccumulator = (192UL << 24); // Start at exactly 270 degrees (Absolute Trough)
  envIncrement = calcPhaseInc(envelopeFreq);

  fadeInterrupts = 0; 
  fadeOutStartInterrupt = 4294967295UL; 

  // Switch to Master Clock
  unsigned long totalInterrupts = (unsigned long)((duration / 1000.0) * actualPWMFreq);
  currentTick = 0; 
  completedCycles = 0;

  if (useChA) { setChA(MidLumi); }  
  if (useChB) { setChB(MidLumi); }  

  PORTC |= (1 << PORTC6);     
  TIMSK0 &= ~_BV(TOIE0); 
  
  setTimer1Callback(sinewaveEnvelopeInterrupt);
  startTimer1Interrupt();

  // ATOMIC READ: Master Clock
  unsigned long safeTick = 0;
  while (true) {
    noInterrupts();
    safeTick = currentTick;
    interrupts();
    
    if (safeTick >= totalInterrupts) break;
    delayMicroseconds(1);  
  }
  
  stopTimer1Interrupt();    
  TIMSK0 |= _BV(TOIE0); 

  PORTD &= ~(1 << PIND4);   
  PORTC &= ~(1 << PORTC6);  
  Serial.print(F("-1\n"));
  Serial.flush();
  
  if (useChA) { setChA(TopLumi / 2); }  
  if (useChB) { setChB(TopLumi / 2); }  
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
      PORTD ^= (1 << PIND4); 
  }

  uint8_t indexA = phaseAccumulatorA >> 24;
  uint8_t indexB = phaseAccumulatorB >> 24;
  uint8_t indexEnv = envAccumulator >> 24;

  uint16_t contrastMultInt = sineWaveTable[indexEnv];  // Localized and optimized math

  uint16_t currentContrastIntA = ((uint32_t)contrastMultInt * contrastMultIntA) >> 8;
  uint16_t currentContrastIntB = ((uint32_t)contrastMultInt * contrastMultIntB) >> 8;

  long tempA = (long)sineWaveTable[indexA] - MidLumi;
  long ocrValA_calc = MidLumi + ((tempA * currentContrastIntA) >> 8);
  if (ocrValA_calc < 0) ocrValA_calc = 0;
  if (ocrValA_calc > 255) ocrValA_calc = 255;
  uint16_t ocrValA = (uint16_t)ocrValA_calc; 

  long tempB = (long)sineWaveTable[indexB] - MidLumi;
  long ocrValB_calc = MidLumi + ((tempB * currentContrastIntB) >> 8);
  if (ocrValB_calc < 0) ocrValB_calc = 0;
  if (ocrValB_calc > 255) ocrValB_calc = 255;
  uint16_t ocrValB = (uint16_t)ocrValB_calc; 

  if (useChA) { setChA(ocrValA); }  
  if (useChB) { setChB(ocrValB); }  
}

///////////////////////////////////  FREQUENCY SWEEP FUNCTIONS (DDS) //////////////////////////////////

void SteppedFrequencySweep(float startFreq, float endFreq, float stepFreq, int cyclesPerFreq, float phaseA, float phaseB, float contrastA, float contrastB) {

  phaseAccumulatorA = (uint32_t)(phaseA * 4294967296.0);
  phaseAccumulatorB = (uint32_t)(phaseB * 4294967296.0);
  contrastMultIntA = (uint16_t)(contrastA * 256.0);
  contrastMultIntB = (uint16_t)(contrastB * 256.0);

  fadeInterrupts = 0; 
  fadeOutStartInterrupt = 4294967295UL; 

  PORTC |= (1 << PORTC6);  
  PORTD |= (1 << PIND4);   

  float currentFreq = startFreq;
  bool sweepingUp = startFreq <= endFreq;
  stepFreq = abs(stepFreq); 

  phaseIncrementA = calcPhaseInc(currentFreq);
  phaseIncrementB = phaseIncrementA;

  setTimer1Callback(sinewaveInterrupt);
  
  TIMSK0 &= ~_BV(TOIE0); 
  startTimer1Interrupt();

  while ((sweepingUp && currentFreq <= endFreq) || (!sweepingUp && currentFreq >= endFreq)) {
    
    uint32_t newInc = calcPhaseInc(currentFreq);
    noInterrupts(); 
    phaseIncrementA = newInc;
    phaseIncrementB = newInc; // Sweeps apply identically to both channels
    interrupts();
    
    // For stepped sweep, we continue to base execution on cycle count
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

  stopTimer1Interrupt();
  TIMSK0 |= _BV(TOIE0); 
  
  PORTD &= ~(1 << PIND4);   
  PORTC &= ~(1 << PORTC6);  
  Serial.print(F("-1\n"));

  if (useChA) { setChA(TopLumi / 2); } 
  if (useChB) { setChB(TopLumi / 2); } 
}

void FrequencySweep(float fmin, float fmax, float sweepFactorPerSec,
                    float phaseA, float phaseB, float contrastA, float contrastB) {

  if (fmin <= 0.0f) fmin = 0.001f;
  if (fmax <= 0.0f) fmax = 0.001f;

  if (fmin > fmax) {
    float temp = fmin;
    fmin = fmax;
    fmax = temp;
  }

  phaseAccumulatorA = (uint32_t)(phaseA * 4294967296.0);
  phaseAccumulatorB = (uint32_t)(phaseB * 4294967296.0);
  contrastMultIntA = (uint16_t)(contrastA * 256.0);
  contrastMultIntB = (uint16_t)(contrastB * 256.0);

  fadeInterrupts = 0; 
  fadeOutStartInterrupt = 4294967295UL;

  unsigned long totalDurationUs = 0;
  bool isSweeping = false;

  // Optimized continuous math logic
  if (sweepFactorPerSec > 0.000001f && fmax > fmin && (fmax / fmin) > 1.000001f) {
    float timeForOneWaySweepSec_calc = log(fmax / fmin) / sweepFactorPerSec;

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

  PORTC |= (1 << PORTC6);  
  PORTD |= (1 << PIND4);   
  
  TIMSK0 &= ~_BV(TOIE0); 
  
  // Master Clock upgrade for Frequency Sweep!
  unsigned long totalInterrupts = (unsigned long)((totalDurationUs / 1000000.0) * actualPWMFreq);
  currentTick = 0;
  completedCycles = 0;
  
  setTimer1Callback(sinewaveInterrupt);
  startTimer1Interrupt();

  unsigned long safeTick = 0;
  
  while (true) {
      noInterrupts();
      safeTick = currentTick;
      interrupts();
      
      if (safeTick >= totalInterrupts) break;
    
    float sweepProgress = (float)safeTick / (float)totalInterrupts;
    
    if (isSweeping) {
      if (sweepProgress < 0.5) {  
        actual_calc_freq = fmin * exp(sweepFactorPerSec * (sweepProgress * totalDurationUs / 1000000.0));
        if (actual_calc_freq > fmax) actual_calc_freq = fmax;
      } else {                    
        actual_calc_freq = fmax * exp(-sweepFactorPerSec * ((sweepProgress - 0.5) * totalDurationUs / 1000000.0));
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

  stopTimer1Interrupt(); 
  TIMSK0 |= _BV(TOIE0); 
     
  PORTD &= ~(1 << PIND4);   
  PORTC &= ~(1 << PORTC6);  
  Serial.print(F("-1\n"));

  if (useChA) { setChA(TopLumi / 2); }  
  if (useChB) { setChB(TopLumi / 2); }  
}

/////////////////////////////////// SOME GENERIC PWM FUNCTIONS ///////////////////////////////////////////

void setChA(uint16_t ocrValue) {
  OCR1A = pgm_read_word_near(currentChALUT + ocrValue);
}

void setChB(uint16_t ocrValue) {
  OCR1B = pgm_read_word_near(currentChBLUT + ocrValue);
}

void setDutyCycle(float dutyCyclePercentage_A, float dutyCyclePercentage_B, long TopLumi) {
  if (dutyCyclePercentage_A < 0.0) dutyCyclePercentage_A = 0.0;
  if (dutyCyclePercentage_A > 100.0) dutyCyclePercentage_A = 100.0;
  uint16_t ocrValueA = (long)((dutyCyclePercentage_A / 100.0) * TopLumi);

  if (dutyCyclePercentage_B < 0.0) { dutyCyclePercentage_B = 0.0; };
  if (dutyCyclePercentage_B > 100.0) { dutyCyclePercentage_B = 100.0; }
  uint16_t ocrValueB = (long)((dutyCyclePercentage_B / 100.0) * TopLumi);

  PORTD ^= (1 << PIND4);  
  setChA(ocrValueA);
  setChB(ocrValueB);
}

void setDutyCycleTime(float dutyCyclePercentage_A, float dutyCyclePercentage_B, long duration, long TopLumi) {
  if (dutyCyclePercentage_A < 0.0) dutyCyclePercentage_A = 0.0;
  if (dutyCyclePercentage_A > 100.0) dutyCyclePercentage_A = 100.0;
  uint16_t ocrValueA = (long)((dutyCyclePercentage_A / 100.0) * TopLumi);

  if (dutyCyclePercentage_B < 0.0) { dutyCyclePercentage_B = 0.0; };
  if (dutyCyclePercentage_B > 100.0) { dutyCyclePercentage_B = 100.0; }
  uint16_t ocrValueB = (long)((dutyCyclePercentage_B / 100.0) * TopLumi);

  long startTime = millis();          
  PORTC |= (1 << PORTC6);             
  PORTD |= (1 << PIND4);              
  if (useChA) { setChA(ocrValueA); }  
  if (useChB) { setChB(ocrValueB); }  
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  
  }
  PORTD &= ~(1 << PIND4);   
  PORTC &= ~(1 << PORTC6);  
  Serial.print(F("-1\n"));
  Serial.flush();
  if (useChA) { setChA(TopLumi / 2); }  
  if (useChB) { setChB(TopLumi / 2); }  
}

void cycleDutyCycles(float stepSize, float waitTime, int nReps, long TopLumi) {
  for (int irep = 0; irep < nReps; irep++) {
    float dutyCycle = 0;
    while (dutyCycle <= 1) {
      Serial.print(dutyCycle);
      Serial.print(F("\n"));
      long ocrValue = (long)(dutyCycle * TopLumi);
      if (useChA) { setChA(ocrValue); }  
      if (useChB) { setChB(ocrValue); }  
      dutyCycle = dutyCycle + stepSize;
      delay(waitTime);
    }
  }
  Serial.print(F("-1\n"));
  if (useChA) { setChA(TopLumi / 2); }  
  if (useChB) { setChB(TopLumi / 2); }  
}

void readAnalogVals() {
  bool keepReading = true;
  const unsigned long interval = 100;  
  unsigned long previousMillis = millis();

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
      if (analogValue0 < 420) {
        currentChALUT = ChA2LUT;
        currentChBLUT = ChA2LUT;
        Serial.print(F("LUT 2 SELECTED\n"));
      } else {
        currentChALUT = ChA1LUT;
        currentChBLUT = ChA1LUT;
        Serial.print(F("LUT 1 SELECTED\n"));
      }
    }
  }
}

////////////////////// TIMER 1 PWM FREQUENCY CONTROL //////////////////////////////

long calculatePrescalerAndTOP(long desiredFrequency, long &prescaler) {
  long TOP = 0;
  long possiblePrescalers[] = { 1, 8, 64, 256, 1024 };

  for (int i = 0; i < 5; i++) {
    long currentPrescaler = possiblePrescalers[i];
    long calculatedTOP = (CLOCK_FREQ / (2 * currentPrescaler * desiredFrequency)) - 1;

    if (calculatedTOP >= 0 && calculatedTOP <= 65535) {
      prescaler = currentPrescaler;
      TOP = calculatedTOP;
      break;  
    }
  }
  return TOP;
}

void configureTimer1(long prescaler, long TOP) {
  TCCR1A = 0;
  TCCR1B = 0;

  TCCR1A |= (1 << COM1A1);  
  TCCR1A |= (1 << COM1B1);  
  TCCR1B |= (1 << WGM13);   
  TCCR1B &= ~(1 << WGM12);  
  TCCR1A &= ~(1 << WGM11);  
  TCCR1A &= ~(1 << WGM10);  

  switch (prescaler) {
    case 1:
      TCCR1B |= (1 << CS10);  
      break;
    case 8:
      TCCR1B |= (1 << CS11);  
      break;
    case 64:
      TCCR1B |= (1 << CS11) | (1 << CS10);  
      break;
    case 256:
      TCCR1B |= (1 << CS12);  
      break;
    case 1024:
      TCCR1B |= (1 << CS12) | (1 << CS10);  
      break;
    default:
      break;
  }

  ICR1 = TOP;
}

////////////////////////////// Timer 1 Overflow Interrupt control /////////////////////////////////

void (*timer1Callback)() = nullptr;  

void setTimer1Callback(void (*callback)()) {
  timer1Callback = callback;
}

void startTimer1Interrupt() {
  TIFR1 |= (1 << TOV1);     
  TIMSK1 |= (1 << TOIE1);   
}

void stopTimer1Interrupt() {
  TIMSK1 &= ~(1 << TOIE1);  
}

ISR(TIMER1_OVF_vect) {
  if (timer1Callback) {
    timer1Callback();  
  }
}