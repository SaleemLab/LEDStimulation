//#include <Firmata.h>

//#define SERIAL_TX_BUFFER_SIZE 2048
//#define SERIAL_RX_BUFFER_SIZE 2048

//#define wheelPin_1 5        // Wheel
//#define wheelPin_2 3        // Wheel
//#define wheelPin_1 2       // Wheel (Grey)
//#define wheelPin_2 5        // Wheel (green)
#define wheelPin_1 20       // Wheel (Grey)
#define wheelPin_2 21       // Wheel (green)
#define cameraPin 7             // digital pin of sync
#define cameraPin2 11             // digital pin of sync
#define syncPin 3             // digital pin of sync

/*
int photodiodePin = 0;    
int syncPin = 0;
*/

volatile long encoderCount = 0;    // Stores the count from the encoder
volatile uint8_t lastEncoded = 0;  // Holds the previous state of the encoder

unsigned int CameraVal;
unsigned int CameraVal2;




unsigned long lastTime = millis();
unsigned long sampleTime = millis();
unsigned int toggle = 0 ;

String msg;
byte Buf[20];



void setup()

{
  pinMode(wheelPin_1, INPUT_PULLUP);   // wheel1
  pinMode(wheelPin_2, INPUT_PULLUP);   // wheel2
  pinMode(cameraPin, INPUT);
  pinMode(cameraPin2, INPUT);
  pinMode(syncPin,OUTPUT);

  digitalWrite(syncPin,LOW);

//  
  Serial.begin(115200);          //  setup serial
  
  
//  Firmata.begin(115200);
  

  delay(500);

     // Read the initial state of the encoder pins
  uint8_t aState = digitalRead(wheelPin_1);
  uint8_t bState = digitalRead(wheelPin_2);
  lastEncoded = (aState << 1) | bState;
  
  // Attach interrupts to both encoder channels
  attachInterrupt(digitalPinToInterrupt(wheelPin_1), updateEncoder, CHANGE);
  attachInterrupt(digitalPinToInterrupt(wheelPin_2), updateEncoder, CHANGE);
}



void loop()

{
  CameraVal = analogRead(cameraPin);
  CameraVal2 = analogRead(cameraPin2);


  if ((millis()-lastTime)>random(20,50))
  {
    toggle = !toggle;
    digitalWrite(syncPin,toggle);
    lastTime  = millis();
  }

    long count;
  //noInterrupts();
  count = encoderCount;
  //interrupts();

  
  //msg = String(wheelVal1) + "\t" + String(wheelVal2) + "\t" + String(CameraVal) + "\t" + String(CameraVal2) + "\t" + String(toggle) + "\t" + String(millis()) + ",";
  msg = String(encoderCount) + "\t" + String(CameraVal) + "\t" + String(CameraVal2) + "\t" + String(toggle) + "\t" + String(millis()) + ",";

  //msg = String(counter) + "\t" + String(wheelVal1) + "\t" + String(wheelVal2) + "\t" + String(Acounter) + "\t" + String(Bcounter) + "\t" + String(millis()) + ",";
  //msg = String(counter) + "\t" + String(wheelVal1) + "\t" + String(wheelVal2) + "\t" + String(CameraVal2) + "\t" + String(toggle) + "\t" + String(millis()) + ",";
//  msg.getBytes(Buf, 20);
  /*
  Serial.print(wheelVal1);
  Serial.print("\t");
  Serial.print(wheelVal2);
  Serial.print("\t");
  Serial.print(CameraVal);
  Serial.print("\t");
  Serial.print(CameraVal2);
  Serial.print("\t");
  Serial.print(toggle);
  Serial.println(",");
  */
//  
Serial.println(msg);
//  Serial.println((char*)Buf);
//Firmata.sendString(Buf);
   
  
}

void updateEncoder() {
  //Serial.println("Do");
  // Read current state of encoder pins
  uint8_t aState = digitalRead(wheelPin_1);
  uint8_t bState = digitalRead(wheelPin_2);
  uint8_t encoded = (aState << 1) | bState;
  
  // Combine previous and current state to determine direction
  int8_t sum = (lastEncoded << 2) | encoded;
  
  // These patterns represent the valid transitions.
  // Depending on the direction of rotation, we update the count:
  if(sum == 0b0001 || sum == 0b0111 || sum == 0b1110 || sum == 0b1000)
    encoderCount++;
  else if(sum == 0b0010 || sum == 0b0100 || sum == 0b1101 || sum == 0b1011)
    encoderCount--;
  
  lastEncoded = encoded;  // Save the current state for the next update
}
