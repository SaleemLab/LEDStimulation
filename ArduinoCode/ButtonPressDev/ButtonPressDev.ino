// Define the pins
const int leftPin = 2;
const int rightPin = 3;

// Variables to track the LEFT button
int leftState = HIGH;         
int lastLeftState = HIGH;     
unsigned long lastDebounceLeft = 0;  

// Variables to track the RIGHT button
int rightState = HIGH;        
int lastRightState = HIGH;    
unsigned long lastDebounceRight = 0; 

// The debounce delay (50ms is good for heavy switches)
unsigned long debounceDelay = 50;    

void setup() {
  Serial.begin(9600);
  
  // Use internal pull-up resistors for both pins
  pinMode(leftPin, INPUT_PULLUP);
  pinMode(rightPin, INPUT_PULLUP);
}

void loop() {
  // 1. Read the current physical state of both pins
  int readingLeft = digitalRead(leftPin);
  int readingRight = digitalRead(rightPin);

  // ==========================================
  // 2. DEBOUNCE LOGIC FOR THE LEFT BUTTON
  // ==========================================
  
  // If the switch changed (due to noise or pressing), reset the timer
  if (readingLeft != lastLeftState) {
    lastDebounceLeft = millis(); 
  }

  // If the pin has been stable for longer than our 50ms delay...
  if ((millis() - lastDebounceLeft) > debounceDelay) {
    
    // ...and if that stable state is different from the last official state:
    if (readingLeft != leftState) {
      leftState = readingLeft;
      
      // If the new stable state is LOW (Pressed)
      if (leftState == LOW) {
        Serial.println("0");
      }
    }
  }

  // ==========================================
  // 3. DEBOUNCE LOGIC FOR THE RIGHT BUTTON
  // ==========================================
  
  // If the switch changed (due to noise or pressing), reset the timer
  if (readingRight != lastRightState) {
    lastDebounceRight = millis(); 
  }

  // If the pin has been stable for longer than our 50ms delay...
  if ((millis() - lastDebounceRight) > debounceDelay) {
    
    // ...and if that stable state is different from the last official state:
    if (readingRight != rightState) {
      rightState = readingRight;
      
      // If the new stable state is LOW (Pressed)
      if (rightState == LOW) {
        Serial.println("1");
      }
    }
  }

  // ==========================================
  // 4. SAVE THE STATES FOR THE NEXT LOOP
  // ==========================================
  lastLeftState = readingLeft;
  lastRightState = readingRight;
}