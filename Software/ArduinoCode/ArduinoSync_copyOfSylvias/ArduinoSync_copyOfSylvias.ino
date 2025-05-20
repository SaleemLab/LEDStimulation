//

#include <Firmata.h>

//#define SERIAL_TX_BUFFER_SIZE 2048
//#define SERIAL_RX_BUFFER_SIZE 2048

//#define wheelPin_1 5        // Wheel
//#define wheelPin_2 3        // Wheel
#define wheelPin_1 2       // Wheel
#define wheelPin_2 5        // Wheel
#define cameraPin 7             // digital pin of sync
#define cameraPin2 11             // digital pin of sync
#define syncPin 3             // digital pin of sync

/*
int photodiodePin = 0;    
int syncPin = 0;
*/


unsigned int wheelVal1;
unsigned int wheelVal2;
unsigned int CameraVal;
unsigned int CameraVal2;
bool lastStateA;
bool lastStateB;
bool changeSigA;
bool changeSigB;
int counter = 0;
int Acounter = 0;
int Bcounter = 0;

unsigned long lastTime = millis();
unsigned long sampleTime = millis();
unsigned int toggle = 0 ;

String msg;
byte Buf[20];



void setup()

{
  pinMode(wheelPin_1, INPUT);   // PhotoDiode
  pinMode(wheelPin_2, INPUT);   // digital pin of sync
  pinMode(cameraPin, INPUT);
  pinMode(cameraPin2, INPUT);
  pinMode(syncPin,OUTPUT);

  digitalWrite(syncPin,LOW);

//  
  Serial.begin(115200);          //  setup serial
  
  
//  Firmata.begin(115200);
  
 lastStateA = digitalRead(wheelPin_1);
 lastStateB = digitalRead(wheelPin_2);

 changeSigA = !lastStateA;
 changeSigB = !lastStateB;
  delay(500);
}



void loop()

{

  wheelVal1 = digitalRead(wheelPin_1);    // read the input pin
  wheelVal2 = digitalRead(wheelPin_2);
  CameraVal = analogRead(cameraPin);
  CameraVal2 = analogRead(cameraPin2);

if (wheelVal1==changeSigA && lastStateA!=wheelVal1)
{
  Acounter++;
  if (Acounter>Bcounter){
    counter++;
    }
   if (Acounter==Bcounter){
    counter--;
  }
  
//  
  }
  

 if (wheelVal2==changeSigB && lastStateB!=wheelVal2)
{
 Bcounter++;
 
} 

if (Bcounter<Acounter-2)
{
  Bcounter=Acounter;
}

if (Acounter<Bcounter-2)
{
  Acounter=Bcounter;
}

  

//  if ((lastStateA!=wheelVal1) && wheelVal1==changeSigA){
//    if (wheelVal2!=changeSigB){
//      counter++;
//    }
//    else
//    {
//      if (wheelVal2==changeSigB){
//      counter--;
//      }
//      }
//  }





  lastStateA = wheelVal1;
  lastStateB = wheelVal2;

  if ((millis()-lastTime)>random(20,50))
  {
    toggle = !toggle;
    digitalWrite(syncPin,toggle);
    lastTime  = millis();
  }

  
  msg = String(wheelVal1) + "\t" + String(wheelVal2) + "\t" + String(CameraVal) + "\t" + String(CameraVal2) + "\t" + String(toggle) + "\t" + String(millis()) + ",";
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
