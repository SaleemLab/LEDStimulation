const int ledPin =  3;// the number of the LED pin
const unsigned long updateTime = 100; // in ms
const unsigned long contrastSwitchTime = 3000; // in ms

// Variables will change:
long randNumber;      // the variable which is supposed to hold the random number
volatile byte contrast = LOW;

// Generally, you should use "unsigned long" for variables that hold time
// The value will quickly become too large for an int to store
unsigned long previousMillis = 0;        // will store last time LED was updated
unsigned long switchPreviousMillis = 0;        // will store last time LED was updated


void setup() {
  // set the digital pin as output:
  pinMode(ledPin, OUTPUT);
  randomSeed(analogRead(0));
  
}

void loop() {


// check if high or low contrast
  unsigned long switchCurrentMillis = millis();
    if (switchCurrentMillis - switchPreviousMillis >= contrastSwitchTime) {
    switchPreviousMillis = switchCurrentMillis;

        // if the LED is off turn it on and vice-versa:
    contrast=!contrast;
    }
  
  unsigned long currentMillis = millis();

  if (currentMillis - previousMillis >= updateTime) {
    // save the last time you blinked the LED
    previousMillis = currentMillis;

  if (contrast==HIGH) 
  {

// set the LED with the ledState of the variable:
     randNumber = random(0, 255);
     analogWrite(ledPin, randNumber);         
  } else
       randNumber = random(64, 191);
     analogWrite(ledPin, randNumber);  
  }
}
