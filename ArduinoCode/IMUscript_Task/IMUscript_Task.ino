// ==========================================
// IMU & Hardware Sync Interrupt Script
// Board: Arduino Leonardo
// Sensor: SparkFun/Adafruit BNO085
// ==========================================

#include <Wire.h>
#include "SparkFun_BNO080_Arduino_Library.h"

// Instantiate the sensor
BNO080 myIMU;

// --- Sync Pin Definitions (Leonardo) ---
const int stimStatePin = 0;   // Channel 1: Stim ON/OFF
const int stimChangePin = 1;  // Channel 2: stimChange toggling

// --- Volatile Variables for ISRs ---
// Declared 'volatile' because they are modified inside an interrupt
volatile unsigned long isr_lastStimOnTime = 0;
volatile unsigned long isr_lastStimOffTime = 0;
volatile unsigned long isr_lastStimChangeTime = 0;

void setup() {
  Serial.begin(115200);
  
  // Wait for Leonardo's USB serial to connect
  while (!Serial) delay(10); 
  
  // Configure the sync pins as inputs
  pinMode(stimStatePin, INPUT); 
  pinMode(stimChangePin, INPUT);

  // Attach the Interrupts mapping D1 and D0 to internal registers
  attachInterrupt(digitalPinToInterrupt(stimStatePin), handleStimState, CHANGE);
  attachInterrupt(digitalPinToInterrupt(stimChangePin), handleStimChange, CHANGE);
  
  Wire.begin();

  Serial.println(F("Initializing BNO085 for Neon Sync..."));

  // Adafruit breakout boards default to I2C address 0x4A
  if (myIMU.begin(0x4A, Wire) == false) {
    Serial.println(F("Failed to find BNO085 chip. Check wiring!"));
    while (1) { delay(10); } // Halt
  }

  Serial.println(F("BNO085 Connected Successfully!"));

  // Enable Linear Acceleration at 10ms (100Hz)
  myIMU.enableLinearAccelerometer(10);
  
  // Enable Game Rotation Vector at 10ms (100Hz)
  myIMU.enableGameRotationVector(10);

  // Print a combined CSV Header including the new sync columns
  Serial.println(F("Time_us,AccelX,AccelY,AccelZ,QuatI,QuatJ,QuatK,QuatReal,StimON_ts,StimOFF_ts,StimChange_ts"));
}

void loop() {
  // dataAvailable() returns true when the sensor has processed new data
  if (myIMU.dataAvailable() == true) {
    
    // 1. Create local variables to hold the safe copies of our timestamps
    unsigned long safe_StimOnTime;
    unsigned long safe_StimOffTime;
    unsigned long safe_StimChangeTime;

    // 2. Briefly pause interrupts to copy the 4-byte variables safely
    // This prevents corrupted data if an interrupt fires mid-copy
    noInterrupts(); 
    safe_StimOnTime = isr_lastStimOnTime;
    safe_StimOffTime = isr_lastStimOffTime;
    safe_StimChangeTime = isr_lastStimChangeTime;
    interrupts(); // Turn interrupts back on immediately

    // 3. Print the main microsecond timestamp
    Serial.print(micros());
    Serial.print(F(","));
    
    // 4. Print Linear Acceleration (m/s^2)
    Serial.print(myIMU.getLinAccelX(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getLinAccelY(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getLinAccelZ(), 4);
    Serial.print(F(","));
    
    // 5. Print Game Rotation Vector (Quaternions)
    Serial.print(myIMU.getQuatI(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getQuatJ(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getQuatK(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getQuatReal(), 4);
    Serial.print(F(","));

    // 6. Print Hardware Sync Timestamps
    Serial.print(safe_StimOnTime);
    Serial.print(F(","));
    Serial.print(safe_StimOffTime);
    Serial.print(F(","));
    Serial.println(safe_StimChangeTime);
  }
}

// ==========================================
//        INTERRUPT SERVICE ROUTINES
// ==========================================

void handleStimState() {
  // Using micros() to match the IMU's timescale
  unsigned long currentTime = micros(); 
  
  if (digitalRead(stimStatePin) == HIGH) {
    isr_lastStimOnTime = currentTime;
  } else {
    isr_lastStimOffTime = currentTime;
  }
}

void handleStimChange() {
  // Record timestamp of the toggle event
  isr_lastStimChangeTime = micros(); 
}