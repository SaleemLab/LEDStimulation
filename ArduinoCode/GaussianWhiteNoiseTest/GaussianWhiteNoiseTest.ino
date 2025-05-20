// --- Function to generate one pseudo-Gaussian value using CLT ---
// (Same function as provided, generates one N(0,1) approx value per call)
float generateGaussianCLT() {
  const int N = 12; // Number of uniform samples to sum
  long sum_random = 0;
  // Expected mean of N * (M-1)/2 for random(M). Here M=1001.
  // Expected mean = 12 * (1000)/2 = 6000. Let's use the calculation for safety.
  long expected_mean_sum = (long)N * (1001 - 1) / 2; // Correct expected mean

  // Variance of U[0, M-1] is (M^2 - 1)/12. For M=1001, approx M^2/12.
  // Variance of Sum = N * Var(U) = N * (M^2 - 1) / 12
  // Std Dev (Sum) = sqrt(N * (M^2 - 1) / 12)
  // For N=12, Std Dev = sqrt(12 * (1001^2 - 1) / 12) = sqrt(1001^2 - 1) approx 1001
  // So the scale factor is approx 1/1001 to get N(0, 1)
  float scale_factor = 1.0 / sqrt((float)N * (pow(1001.0, 2) - 1.0) / 12.0);
  // Pre-calculated simplification for N=12 is roughly 1.0 / 1000.999... ~= 1.0/1001
  // float scale_factor = 1.0 / 1001.0; // Using the simplified scaling

  for (int i = 0; i < N; i++) {
    sum_random += random(1001); // Generate uniform integer in [0, 1000]
  }

  // Calculate deviation from the mean and scale to approximate N(0, 1)
  float gaussian_approx = ((float)sum_random - (float)expected_mean_sum) * scale_factor;
  return gaussian_approx;
}

// --- Global variables ---
volatile float latest_gaussian_value_1 = 0.0; // Stores the latest value for sequence 1
volatile float latest_gaussian_value_2 = 0.0; // Stores the latest value for sequence 2
volatile bool needs_new_values = true;      // Flag to signal when to calculate new values
unsigned long last_output_time = 0;
const unsigned long output_interval = 8; // Output every 8 milliseconds (~125 Hz)

void setup() {
  Serial.begin(115200); // Ensure your Serial Monitor matches this
  while (!Serial); // Wait for Serial port to connect (especially needed for Leonardo)

  // Seed the random number generator ONCE in setup
  // Use an unconnected analog pin for a somewhat unpredictable seed
  randomSeed(analogRead(A0));

  pinMode(LED_BUILTIN, OUTPUT); // Optional: for visualizing timing
  Serial.println("Setup complete. Starting generation...");
}

void loop() {
  // --- Check if new values need to be calculated ---
  if (needs_new_values) {
    digitalWrite(LED_BUILTIN, HIGH); // Optional: Indicate calculation start
    latest_gaussian_value_1 = generateGaussianCLT(); // Generate first value
    latest_gaussian_value_2 = generateGaussianCLT(); // Generate second value
    needs_new_values = false; // Mark that we have fresh values
    digitalWrite(LED_BUILTIN, LOW); // Optional: Indicate calculation end
  }

  // --- Check if it's time to output/use the latest values ---
  unsigned long current_time = millis();
  if (current_time - last_output_time >= output_interval) {
    last_output_time = current_time;

    // --- Use the latest generated gaussian values here ---
    // Example: print them comma-separated
    Serial.print(latest_gaussian_value_1, 6); // Print with 6 decimal places
    Serial.print(",");
    Serial.println(latest_gaussian_value_2, 6); // Print with 6 decimal places

    // ----------------------------------------------------

    needs_new_values = true; // Signal that new values should be calculated next loop cycle

    // Optional: Blink LED to confirm output timing (can be fast!)
    // digitalWrite(LED_BUILTIN,!digitalRead(LED_BUILTIN));
  }

  // The loop continues. If needs_new_values is false, it skips the
  // calculation block until after the next output occurs.
  // Other non-blocking code could potentially run here if needed.
}