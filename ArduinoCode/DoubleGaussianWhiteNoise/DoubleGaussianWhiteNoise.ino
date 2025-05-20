#include <math.h> // Include for sqrt and pow if using the more precise scale_factor

// --- Configuration ---
const float TARGET_MEAN = 128.0;      // Set your desired mean here
const float TARGET_STD_DEV = 50.0;   // Set your desired standard deviation here
const unsigned long output_interval = 8; // Output every 8 milliseconds (~125 Hz)
const int CLT_N = 12; // Number of uniform samples to sum for CLT

// --- Function to generate one pseudo-Gaussian value using CLT ---
// Returns an approximation of a N(0, 1) value (standard normal distribution)
float generateGaussianCLT() {
  const int N = CLT_N;
  long sum_random = 0;
  // Expected mean of N * (M-1)/2 for random(M). Here M=1001.
  long expected_mean_sum = (long)N * (1001 - 1) / 2; // Expected mean of the sum

  // Variance of U[0, M-1] is (M^2 - 1)/12. For M=1001.
  // Variance of Sum = N * Var(U) = N * (M^2 - 1) / 12
  // Std Dev (Sum) = sqrt(N * (M^2 - 1) / 12)
  // For N=12, Std Dev approx 1001.
  // Scale factor is 1 / Std Dev(Sum) to normalize the sum to N(0, 1)
  // More precise calculation:
  float scale_factor = 1.0 / sqrt((float)N * (pow(1001.0, 2) - 1.0) / 12.0);
  // Simpler approximation often used for N=12:
  // float scale_factor = 1.0 / 1001.0;

  for (int i = 0; i < N; i++) {
    sum_random += random(1001); // Generate uniform integer in [0, 1000]
  }

  // Calculate deviation from the mean and scale to approximate N(0, 1)
  float gaussian_approx_zero_mean_unit_variance =
      ((float)sum_random - (float)expected_mean_sum) * scale_factor;

  return gaussian_approx_zero_mean_unit_variance;
}

// --- Global variables ---
volatile float latest_gaussian_value_1 = 0.0; // Stores the latest scaled value for sequence 1
volatile float latest_gaussian_value_2 = 0.0; // Stores the latest scaled value for sequence 2
volatile bool needs_new_values = true;      // Flag to signal when to calculate new values
unsigned long last_output_time = 0;


void setup() {
  Serial.begin(115200); // Ensure your Serial Monitor matches this
  while (!Serial); // Wait for Serial port to connect (especially needed for Leonardo)

  // Seed the random number generator ONCE in setup
  randomSeed(analogRead(A0)); // Use an unconnected analog pin for a somewhat unpredictable seed

  pinMode(9, OUTPUT); // Optional: for visualizing timing
  pinMode(10, OUTPUT); // Optional: for visualizing timing

  Serial.println("Setup complete.");
  Serial.print("Target Mean: "); Serial.println(TARGET_MEAN);
  Serial.print("Target Std Dev: "); Serial.println(TARGET_STD_DEV);
  Serial.println("Starting generation (Seq1, Seq2):");
}

void loop() {
  // --- Check if new values need to be calculated ---
  if (needs_new_values) {
    // digitalWrite(LED_BUILTIN, HIGH); // Optional: Indicate calculation start

    // Generate first standard normal value Z ~ N(0, 1)
    float z1 = generateGaussianCLT();
    // Scale it to X ~ N(mean, std_dev^2) using X = mean + std_dev * Z
    latest_gaussian_value_1 = TARGET_MEAN + TARGET_STD_DEV * z1;

    // Generate second standard normal value Z ~ N(0, 1)
    float z2 = generateGaussianCLT();
     // Scale it to X ~ N(mean, std_dev^2) using X = mean + std_dev * Z
    latest_gaussian_value_2 = TARGET_MEAN + TARGET_STD_DEV * z2;

    needs_new_values = false; // Mark that we have fresh values
    // digitalWrite(LED_BUILTIN, LOW); // Optional: Indicate calculation end
  }

  // --- Check if it's time to output/use the latest values ---
  unsigned long current_time = millis();
  if (current_time - last_output_time >= output_interval) {
    last_output_time = current_time;

    // --- Use the latest generated scaled gaussian values here ---
    // Example: print them comma-
    //Serial.print(latest_gaussian_value_1, 6); // Print with 6 decimal places
//Serial.print(",");
   // Serial.print(latest_gaussian_value_2, 6); // Print with 6 decimal places
//Serial.print(",");
//Serial.println(current_time);

    // ----------------------------------------------------

    analogWrite(9, latest_gaussian_value_1);
    analogWrite(10, latest_gaussian_value_2);


    needs_new_values = true; // Signal that new values should be calculated next loop cycle
  }
   // Other non-blocking code can run here
}