#define CLOCK_FREQ 16000000  // Arduino Leonardo clock frequency (16 MHz)
#define TABLE_SIZE 256       // Number of samples in the wavetable

long prescaler;
long TOP;      // set by the clock
long TopLumi;  // use to limit max luminance

// serial
const byte numChars = 30;
char receivedChars[numChars];  // an array to store the received data
bool newData = false;

// stimulus selection char
String FirstChar;

// Array to hold the wavetable
int sineWaveTable[TABLE_SIZE];
// wavetable step size
volatile int stepSize = 1;
// contrast-envelope counter and contrast multiplier
volatile int nEnvCounts = 1;
volatile float contrastMult = 1;

// white noise update time and PWM value
volatile long updateTime;
volatile long randNumber = 0;

// flicker state
volatile bool toggleState = false;  // Flag to track the current state

// array for frozen white noise values
volatile int frozenWhiteNoiseTable[375];  // max number of values
volatile int frozenWhiteNoiseTableSize;



void setup() {
  pinMode(9, OUTPUT);      // Pin 9 controlled by Timer1 (Channel A)
  pinMode(4, OUTPUT);      // Pin 4 indicator pin, e.g. for sinewave cycles
  pinMode(5, OUTPUT);      // Pin 5 stim ON or OFF pin 

  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW


  Serial.begin(115200);

  // Set the desired PWM frequency (in Hz)
  long desiredPWMFrequency = 7680;  // Example: 7680 Hz

  // Constrain the desired PWM frequency to be a multiple of TABLE_SIZE
  // (probably not important since precision seems low at high frequenies)
  desiredPWMFrequency = constrainFrequency(desiredPWMFrequency);

  // Calculate the prescaler and TOP value for the constrained frequency
  TOP = calculatePrescalerAndTOP(desiredPWMFrequency, prescaler);

  // Apply the prescaler and TOP value to Timer1
  configureTimer1(prescaler, TOP);

  TopLumi = TOP;  // set default TopLumi (max duty cycle) as the actual TOP (i.e. 100% for now)

  // Generate the sine wave LUT based on the Timer1 config and TopLumi
  generateSineWaveTable(TopLumi);

  // initialise random number to 50% duty cyle for white noise stimuli
  randNumber = TopLumi / 2;

  // Set pin 9 to 50% duty cycle as default
  OCR1A = TopLumi / 2;

  //delay(5000);
  //whiteNoise(10000, 10);
  //outputSinewave(10,10000);

  //ToggleLEDTest(8000);
}

// loop runs checking for new serial input
void loop() {
  GetSerialInput();
  if (newData) {
    newData = false;
    ActionSerial();
  }
}



//////////////////////////////////////// HANDLE SERIAL INPUT //////////////////////////////////////
void GetSerialInput() {  // part of code taken from http://forum.arduino.cc/index.php?topic=396450.0
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
    } else {                      // serial message finished
      receivedChars[ndx] = '\0';  // terminate the string
      ndx = 0;
      newData = true;
    }
  }
}

void ActionSerial() {  // Actions serial data by choosing appropriate stimulation
  Serial.print("rc: ");
  Serial.println(receivedChars);
  char delimiters[] = ",";
  char *token;
  uint8_t idx = 0;
#define MAX_VALS 5  // max required? freq, duration, contrast, carrier freq?
  char *serialVals[MAX_VALS];
  token = strtok(receivedChars, ",");


  while (token != NULL) {
    //Serial.println( token );
    if (idx < MAX_VALS)
      serialVals[idx++] = token;
    token = strtok(NULL, ",");
  }

  FirstChar = serialVals[0];

  if (FirstChar == "s")  // sinusoidal flicker
  {
    long stimulusDuration = atof(serialVals[1]);
    float frequency = atof(serialVals[2]);
    Serial.println("Stim: Sinusoidal dimming");
    Serial.print("Stim duration: ");
    Serial.println(stimulusDuration);
    Serial.print("Frequency: ");
    Serial.println(frequency);

    outputSinewave(frequency, stimulusDuration);
  } else if (FirstChar == "wn")  // white noise
  {
    long stimulusDuration = atof(serialVals[1]);
    long updateTime = atof(serialVals[2]);
    Serial.println("Stim: White noise");
    Serial.print("Stim duration: ");
    Serial.println(stimulusDuration);
    Serial.print("Update time: ");
    Serial.println(updateTime);

    whiteNoise(updateTime, stimulusDuration);
  } else if (FirstChar == "fwn")  // frozen white noise
  {
    long stimulusDuration = atof(serialVals[1]);
    int updateTime = atof(serialVals[2]);
    int nReps = atof(serialVals[3]);
    Serial.println("Stim: White noise");
    Serial.print("Stim duration: ");
    Serial.println(stimulusDuration);
    Serial.print("Update time: ");
    Serial.println(updateTime);

    frozenWhiteNoise(updateTime, stimulusDuration,nReps);
  } else if (FirstChar == "se")  // sinusoidal flicker with contrast envelope
  {
    long stimulusDuration = atof(serialVals[1]);
    float frequency = atof(serialVals[2]);
    float envFrequency = atof(serialVals[3]);
    Serial.println("Stim: Sinusoidal env");
    Serial.print("Stim duration: ");
    Serial.println(stimulusDuration);
    Serial.print("Frequency: ");
    Serial.println(frequency);
    Serial.print("Envelope freq: ");
    Serial.println(envFrequency);
    SineContrastConv(stimulusDuration, frequency, envFrequency);
  } else if (FirstChar == "f")  // Square wave flicker
  {
    long stimulusDuration = atof(serialVals[1]);
    float FlickerFreq = atof(serialVals[2]);
    Serial.println("Stim: Square Wave flicker");
    Serial.print("Stim duration: ");
    Serial.println(stimulusDuration);
    Serial.print("Frequency: ");
    Serial.println(FlickerFreq);
    FlickerLED(FlickerFreq, stimulusDuration);
  } else if (FirstChar == "st")  // Set TopLuminance
  {
    float TopMultiplier = atof(serialVals[1]);

    SetTopLumi(TopMultiplier);

  } else if (FirstChar == "sd")  // Set duty cycle
  {
    float DutyCyle = atof(serialVals[1]);

    setDutyCycle(DutyCyle, TopLumi);

  } else  // not valid stimulus code
  {
    Serial.print(FirstChar);
    Serial.println(" is an invalid stimulus code - make sure you are using carriage return line ending");
  }
  //memset('\0', receivedChars, sizeof(receivedChars));
  memset(receivedChars, '\0', sizeof(receivedChars));
  //Serial.print("rc: ");
  //Serial.println(receivedChars);
}


///////////////////////////////////// SINEWAVE FLICKER  //////////////////////////////////////
// Function to generate a sine wave table
void generateSineWaveTable(long TOP) {
  for (int i = 0; i < TABLE_SIZE; i++) {
    // Calculate the sine wave value (scaled between 0 and TOP)
    float angle = (2.0 * PI * i) / TABLE_SIZE;                     // Angle in radians
    sineWaveTable[i] = (int)((sin(angle) + 1.0) * (TOP / 2.0));  // Scale to 0-TOP
  }
}

void outputSinewave(float sinewaveFrequency, long duration) {

  int tableIndex = 0;  // Start at the beginning of the sine wave table

  // first do some calculations to find the update interval and step size for Timer3 interrupts
  // Calculate the PWM cycle time in microseconds
  float pwmCycleTime = (2.0 * TOP) / (float)(CLOCK_FREQ / prescaler);  // Time per PWM cycle in seconds
  pwmCycleTime *= 1e6;                                                 // Convert seconds to microseconds

  // Calculate the base update interval for the sinewave frequency
  float baseUpdateInterval = 1.0 / (sinewaveFrequency * TABLE_SIZE);  // Time per table update in seconds
  baseUpdateInterval *= 1e6;                                          // Convert to microseconds

  // Find the smallest step size that is a factor of TABLE_SIZE
  stepSize = TABLE_SIZE;  // Start with the maximum possible step size
  for (int i = 1; i <= TABLE_SIZE; i++) {
    if ((TABLE_SIZE % i == 0) && (baseUpdateInterval * i >= pwmCycleTime)) {
      stepSize = i;
      break;  // Stop at the first valid (smallest) step size
    }
  }

  // Recalculate the effective update interval based on the step size
  float updateInterval = baseUpdateInterval * stepSize;
  float updateFrequency = 1e6 / updateInterval;  // update frequency for timer3 interrupt

  // Configure timer3 interrupt to updateFrequency
  configureTimer3Interrupt(updateFrequency);


  long startTime = millis();  // Record the start time
  PORTC |= (1 << PORTC6);
  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(sinewaveInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // finish playing sinewave
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");

  OCR1A = TopLumi / 2;
}

// sinewave interrupt function
void sinewaveInterrupt() {
  static int tableIndex = 0;  // Start at the beginning of the sine wave table
  // Update PWM duty cycle with the next sine wave value
  OCR1A = sineWaveTable[tableIndex];
  if (tableIndex == 0) {
    PORTD ^= (1 << PIND4);  // Toggle Pin 4 if tableIndex is 0
  }

  // Update the table index (wrap around if necessary)
  tableIndex = (tableIndex + stepSize) % TABLE_SIZE;
}




/////////////////////////////////// SINE WAVE FLICKER WITH CONTRAST ENVELOPE //////////////////////////
void SineContrastConv(float duration, float sinewaveFrequency, float envelopeFreq) {

  int tableIndex = 0;  // Start at the beginning of the sine wave table

  // first do some calculations to find the update interval and step size for Timer3 interrupts
  // Calculate the PWM cycle time in microseconds
  float pwmCycleTime = (2.0 * TOP) / (float)(CLOCK_FREQ / prescaler);  // Time per PWM cycle in seconds
  pwmCycleTime *= 1e6;                                                 // Convert seconds to microseconds

  // Calculate the base update interval for the sinewave frequency
  float baseUpdateInterval = 1.0 / (sinewaveFrequency * TABLE_SIZE);  // Time per table update in seconds
  baseUpdateInterval *= 1e6;                                          // Convert to microseconds

  // Find the smallest step size that is a factor of TABLE_SIZE
  stepSize = TABLE_SIZE;  // Start with the maximum possible step size
  for (int i = 1; i <= TABLE_SIZE; i++) {
    if ((TABLE_SIZE % i == 0) && (baseUpdateInterval * i >= pwmCycleTime)) {
      stepSize = i;
      break;  // Stop at the first valid (smallest) step size
    }
  }

  // Recalculate the effective update interval based on the step size
  float updateInterval = baseUpdateInterval * stepSize;
  float updateFrequency = 1e6 / updateInterval;  // update frequency for timer3 interrupt

  // Configure timer3 interrupt to updateFrequency
  configureTimer3Interrupt(updateFrequency);


  // now get interval for updating envelope
  float baseEnvUpdateInterval = 1.0 / (envelopeFreq * TABLE_SIZE);  // Time per table update in seconds
  baseEnvUpdateInterval *= 1e6;                                     // Convert to microseconds

  nEnvCounts = baseEnvUpdateInterval / updateInterval;
  //Serial.print("nc: ");
  //Serial.println(nEnvCounts);

  long startTime = millis();  // Record the start time
  // set timer3 interrupt callback function to play the sinewave
  PORTC |= (1 << PORTC6);  // Stim on pin 5
  setTimer3Callback(sinewaveEnvelopeInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }
  stopTimer3Interrupt();  // finish playing sinewave
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  OCR1A = TopLumi / 2;
}



// sinewave contrast envelope interrupt function
void sinewaveEnvelopeInterrupt() {
  static int tableIndex = 0;  // Start at the beginning of the sine wave table
  static int tableEnvIndex = 0;
  static int envCount = 0;  // contrast envelope counter

  // Update PWM duty cycle with the next sine wave value
  OCR1A = TopLumi / 2 + (sineWaveTable[tableIndex] - TopLumi / 2) * contrastMult;
  if (tableIndex == 0) {
    PORTD ^= (1 << PIND4);  // Toggle Pin 4 if tableIndex is 0
  }

  // Update the table index (wrap around if necessary) for the sinewave carrier
  tableIndex = (tableIndex + stepSize) % TABLE_SIZE;

  // update counter for contrast envelope
  envCount = envCount + 1;
  if (envCount > nEnvCounts - 1)  // if time to update contrast envelope
  {
    envCount = 0;                                                  // reset to 1
    tableEnvIndex = (tableEnvIndex + 1) % TABLE_SIZE;              // get the next contrast value index
    contrastMult = sineWaveTable[tableEnvIndex] / float(TopLumi);  // use index to get the normalised contrast value
    //Serial.println(contrastMult);
  }
}


/////////////////////////////////// WHITE NOISE PWM FUNCTIONS //////////////////////////////

void whiteNoise(long updateTime, long duration) {

  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);


  Serial.print("TOP: ");
  Serial.println(TopLumi);

  long startTime = millis();  // Record the start time
  PORTC |= (1 << PORTC6);  // Stim on pin 5

  // set timer3 interrupt callback function to play the sinewave
  setTimer3Callback(whiteNoiseInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop white noise
  delay(updateTime);
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  OCR1A = TopLumi / 2;
}

// whitenoise interrupt function
void whiteNoiseInterrupt() {
  // Update PWM duty cycle with the next sine wave value
  OCR1A = randNumber;
  PIND = (1 << PIND4);  // alternate PIN 4 value indicator pin
  Serial.println(randNumber);
  randNumber = random(0, TopLumi);  // get a new random number ready
}


/////////////////////////////////// FROZEN WHITE NOISE PWM FUNCTIONS //////////////////////////////

void frozenWhiteNoise(int updateTime, long duration, long nReps) {

  float updateFrequency = 1e3 / updateTime;
  configureTimer3Interrupt(updateFrequency);

  frozenWhiteNoiseTableSize = duration / updateTime;
  long totalDuration = duration * nReps;
  Serial.print("LD: ");
  Serial.println(totalDuration);

  for (int i = 0; i < frozenWhiteNoiseTableSize; i++) {
    // Calculate the sine wave value (scaled between 0 and TOP)
    frozenWhiteNoiseTable[i] = (int)random(0, TopLumi);  // Scale to 0-TOP
  }

  Serial.print("TOP: ");
  Serial.println(TopLumi);
  
  long startTime = millis();  // Record the start time
  // set timer3 interrupt callback function to play the sinewave
  PORTC |= (1 << PORTC6);  // Stim on pin 5
  setTimer3Callback(frozenWhiteNoiseInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < totalDuration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop white noise
  delay(updateTime);
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  OCR1A = TopLumi / 2;
}

// whitenoise interrupt function
void frozenWhiteNoiseInterrupt() {

  static int tableIndex = 0;  // Start at the beginning of the sine wave table
  // Update PWM duty cycle with the frozen white noise value
  OCR1A = frozenWhiteNoiseTable[tableIndex];
  if (tableIndex == 0) {
    PORTD ^= (1 << PIND4);  // Toggle Pin 4 if tableIndex is 0
  }
  //Serial.print("ti: ");
  //Serial.println(tableIndex);
  Serial.println(frozenWhiteNoiseTable[tableIndex]);
  // Update the table index (wrap around at actual white noise table size)
  tableIndex = (tableIndex + 1) % frozenWhiteNoiseTableSize;  //
}

////////////////////////////////////////// SQUARE WAVE FLICKER STIMULUS ////////////////////////////////
void FlickerLED(float flickerFreq, long duration) {
  const long interval = 1e6 / (flickerFreq * 2);  // interval at which to blink (microseconds)
  float updateFrequency = 1e6 / interval;         // update frequency for timer3 interrupt

  configureTimer3Interrupt(updateFrequency);

  long startTime = millis();  // Record the start time
  PORTC |= (1 << PORTC6);  // Stim on pin 5
  setTimer3Callback(SquareWaveFlickerInterrupt);

  // Loop until the specified duration has elapsed
  while (millis() - startTime < duration) {
    delayMicroseconds(1);  //wait for time to end
  }

  stopTimer3Interrupt();  // stop flicker
  PORTD &= ~(1 << PIND4);  // Ensure Pin 4 is set to LOW by changing register directly
  PORTC &= ~(1 << PORTC6); // Ensure Pin 5 is set to LOW
  Serial.println("-1");
  OCR1A = TopLumi / 2;
}

// square wave flicker interrupt function
void SquareWaveFlickerInterrupt() {
  OCR1A = (OCR1A == TopLumi) ? 0 : TopLumi;  // Toggle between 0 and TOP
  PORTD ^= (1 << PIND4);                     // set indicator pin similarly
}


/////////////////////////////////// SOME GENERIC PWM FUNCTIONS ///////////////////////////////////////////
// function to artifically lower the max PWM duty cycle. (i.e. TopMultiplier=0.5 means max duty cycle of 50%)
// other functions will work as normal but scale to this TOP value
void SetTopLumi(float TopMultiplier) {

  if (TopMultiplier > 1) {
    TopMultiplier = 1;
  } else if (TopMultiplier <= 0) {
    TopMultiplier = 1;
  }

  // get new TOP value to use
  TopLumi = float(TOP) * TopMultiplier;

  // Generate the sine wave LUT based on the Timer1 config
  generateSineWaveTable(TopLumi);

  // initialise random number to 50% duty cyle
  randNumber = TopLumi / 2;

  // Set pin 9 to 50% duty cycle as default
  OCR1A = TopLumi / 2;
}


// set the duty cycle manually until a new value is requested
void setDutyCycle(float dutyCyclePercentage, long TopLumi) {
  // Constrain the duty cycle percentage between 0% and 100%
  if (dutyCyclePercentage < 0.0) dutyCyclePercentage = 0.0;
  if (dutyCyclePercentage > 100.0) dutyCyclePercentage = 100.0;

  // Calculate the OCR1A value based on the duty cycle and TOP
  long ocrValue = (long)((dutyCyclePercentage / 100.0) * TopLumi);

  // Set OCR1A to control the duty cycle
  OCR1A = ocrValue;
}

////////////////////// TIMER 1 PWM FREQUENCY CONTROL //////////////////////////////

// Function to calculate the required prescaler and TOP value
long calculatePrescalerAndTOP(long desiredFrequency, long &prescaler) {
  long TOP = 0;

  // Possible prescaler values: 1, 8, 64, 256, 1024
  long possiblePrescalers[] = { 1, 8, 64, 256, 1024 };

  // Try each prescaler and calculate the corresponding TOP
  for (int i = 0; i < 5; i++) {
    long currentPrescaler = possiblePrescalers[i];

    // Calculate the TOP value
    long calculatedTOP = (CLOCK_FREQ / (2 * currentPrescaler * desiredFrequency)) - 1;

    // Check if the calculated TOP value is within the 16-bit range (0 to 65535)
    if (calculatedTOP >= 0 && calculatedTOP <= 65535) {
      prescaler = currentPrescaler;
      TOP = calculatedTOP;
      break;  // Stop after finding the first valid prescaler and TOP
    }
  }

  return TOP;
}

// Function to configure Timer1 with the calculated prescaler and TOP value
void configureTimer1(long prescaler, long TOP) {
  TCCR1A = 0;
  TCCR1B = 0;

  // Set Timer1 in 16-bit Phase Correct PWM mode
  TCCR1A |= (1 << COM1A1);  // Enable PWM on pin 9 (Channel A)
  TCCR1B |= (1 << WGM13);   // Set WGM13 bit
  TCCR1B &= ~(1 << WGM12);  // Clear WGM12 bit
  TCCR1A &= ~(1 << WGM11);  // Clear WGM11 bit
  TCCR1A &= ~(1 << WGM10);  // Clear WGM10 bit

  // Set the prescaler
  switch (prescaler) {
    case 1:
      TCCR1B |= (1 << CS10);  // Prescaler = 1
      break;
    case 8:
      TCCR1B |= (1 << CS11);  // Prescaler = 8
      break;
    case 64:
      TCCR1B |= (1 << CS11) | (1 << CS10);  // Prescaler = 64
      break;
    case 256:
      TCCR1B |= (1 << CS12);  // Prescaler = 256
      break;
    case 1024:
      TCCR1B |= (1 << CS12) | (1 << CS10);  // Prescaler = 1024
      break;
    default:
      break;
  }

  // Set the TOP value
  ICR1 = TOP;
}

// Function to constrain the frequency to a multiple of TABLE_SIZE
long constrainFrequency(long frequency) {
  return frequency - (frequency % TABLE_SIZE);
}


////////////////////////////// Timer3 interrupt control /////////////////////////////////

// Define a function pointer for the interrupt handler
void (*timer3Callback)() = nullptr;  // Initialize to null

// Function to set the callback for Timer3
void setTimer3Callback(void (*callback)()) {
  timer3Callback = callback;
}

// ISR for Timer3 Compare Match A - run the function each interrupt
ISR(TIMER3_COMPA_vect) {
  if (timer3Callback) {
    timer3Callback();  // Call the assigned callback function
  }
}

// Set timer3 frequency and configure for interrupts using CTC mode
void configureTimer3Interrupt(float frequency) {

  long t = micros();
  long prescaler = 0;
  long compareValue = 0;

  // Prescaler options: 1, 8, 64, 256, 1024
  long prescalerOptions[] = { 1, 8, 64, 256, 1024 };
  int prescalerIndex = 0;

  // Iterate through prescaler options to find a valid one
  for (prescalerIndex = 0; prescalerIndex < 5; prescalerIndex++) {
    prescaler = prescalerOptions[prescalerIndex];
    compareValue = (CLOCK_FREQ / (prescaler * frequency)) - 1;

    // Check if compareValue is within the valid 16-bit range
    if (compareValue >= 0 && compareValue <= 65535) {
      break;  // Found a valid prescaler and compareValue
    }
  }

  // If no valid prescaler is found, set to maximum possible values
  if (compareValue < 0 || compareValue > 65535) {
    Serial.println("Unable to configure timer for requested frequency. Adjusting to closest possible.");
    prescaler = 1024;
    compareValue = 65535;
  }

  if (compareValue < 0 || compareValue > 65535) {
    Serial.println("compare value bug!");
    compareValue = 65535;  // Ensure it fits in 16 bits
  }

  // Set Timer3 to CTC mode
  TCCR3A = 0;             // Normal operation
  TCCR3B = (1 << WGM32);  // CTC mode (clear on compare match)

  // Set the prescaler
  TCCR3B &= ~(1 << CS32 | 1 << CS31 | 1 << CS30);  // Clear prescaler bits
  switch (prescaler) {
    case 1: TCCR3B |= (1 << CS30); break;
    case 8: TCCR3B |= (1 << CS31); break;
    case 64: TCCR3B |= (1 << CS31) | (1 << CS30); break;
    case 256: TCCR3B |= (1 << CS32); break;
    case 1024: TCCR3B |= (1 << CS32) | (1 << CS30); break;
  }

  // Set the compare match value
  OCR3A = compareValue;

  // Clear the Timer3 Compare Match A interrupt flag
  TIFR3 |= (1 << OCF3A);  // Writing a 1 clears the flag

  // Enable Timer3 Compare Match A interrupt
  TIMSK3 |= (1 << OCIE3A);

  // Enable global interrupts
  sei();

  //t = micros() - t;
  //Serial.print("micros: ");
  //Serial.println(t);
}

// stop the interrupt function
void stopTimer3Interrupt() {
  // Disable the Timer3 interrupt
  TIMSK3 &= ~(1 << OCIE3A);  // Disable Timer3 Compare Match A interrupt
}