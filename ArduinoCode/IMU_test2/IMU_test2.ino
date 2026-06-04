#include <Wire.h>
#include "SparkFun_BNO080_Arduino_Library.h"

// Instantiate the sensor
BNO080 myIMU;

void setup() {
  Serial.begin(115200);
  
  // Wait for Leonardo's USB serial to connect
  while (!Serial) delay(10); 
  
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

  // Print a combined CSV Header
  Serial.println(F("Time_us,AccelX,AccelY,AccelZ,QuatI,QuatJ,QuatK,QuatReal"));
}

void loop() {
  // dataAvailable() returns true when the sensor has processed new data
  if (myIMU.dataAvailable() == true) {
    
    // 1. Print the microsecond timestamp
    Serial.print(micros());
    Serial.print(F(","));
    
    // 2. Print Linear Acceleration (m/s^2)
    Serial.print(myIMU.getLinAccelX(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getLinAccelY(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getLinAccelZ(), 4);
    Serial.print(F(","));
    
    // 3. Print Game Rotation Vector (Quaternions)
    Serial.print(myIMU.getQuatI(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getQuatJ(), 4);
    Serial.print(F(","));
    Serial.print(myIMU.getQuatK(), 4);
    Serial.print(F(","));
    Serial.println(myIMU.getQuatReal(), 4);
  }
}