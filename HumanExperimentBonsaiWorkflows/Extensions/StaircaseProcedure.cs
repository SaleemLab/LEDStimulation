using System;
using System.Linq; // Added for .Average()
using System.Collections.Generic; // Added for List
using System.Reactive.Linq;
using System.ComponentModel;
using Bonsai;

// We define a custom object to hold the Intensity, Reversal count, and Average
public class StaircaseState
{
    public double Intensity { get; set; }
    public int Reversals { get; set; }
    public double ReversalAverage { get; set; } // NEW: Exposes the rolling average
}

[Description("Adaptive Staircase tracking intensity and reversals with dynamic step sizes and reversal averaging.")]
public class StaircaseProcedure : Transform<bool, StaircaseState>
{
    public StaircaseProcedure()
    {
        InitialIntensity = 10.0;
        StepSizes = new double[] { 2.0, 1.0, 0.5 }; 
        MinIntensity = 0.0;
        MaxIntensity = 100.0;
        CorrectRule = 2;   
        IncorrectRule = 1; 
        MovingAverageWindow = 4; // Default to averaging the last 4 reversals
    }

    [Description("The initial stimulus intensity before any responses are made.")]
    public double InitialIntensity { get; set; }

    [Description("Array of step sizes mapped to reversals (e.g., 0 reversals uses the 1st number, 1 reversal uses the 2nd).")]
    public double[] StepSizes { get; set; } 

    [Description("The number of recent reversals to include in the average calculation.")]
    public int MovingAverageWindow { get; set; } // Property to set N

    [Description("Minimum allowed intensity.")]
    public double MinIntensity { get; set; }

    [Description("Maximum allowed intensity.")]
    public double MaxIntensity { get; set; }

    [Description("Number of consecutive correct responses required to INCREASE intensity.")]
    public int CorrectRule { get; set; }

    [Description("Number of consecutive incorrect responses required to DECREASE intensity.")]
    public int IncorrectRule { get; set; }

    // Helper enum to track the direction the staircase is moving
    private enum Direction { None, Up, Down }

    public override IObservable<StaircaseState> Process(IObservable<bool> source)
    {
        return Observable.Defer(() => {
            
            // Internal state tracking
            double currentIntensity = InitialIntensity;
            int consecutiveCorrect = 0;
            int consecutiveIncorrect = 0;
            int reversals = 0;
            Direction lastDirection = Direction.None;
            
            // A list to hold our recent reversal intensities (peaks and valleys)
            List<double> reversalPoints = new List<double>();

            // source is our incoming stream of booleans (true = correct, false = incorrect)
            return source.Select(isCorrect => {
                bool changed = false;
                Direction currentDirection = Direction.None;
                
                // Capture the intensity BEFORE we step (this is the true peak/valley)
                double preStepIntensity = currentIntensity;

                // --- Determine current step size from the array based on reversals ---
                int safeIndex = (StepSizes != null && StepSizes.Length > 0) 
                                ? Math.Min(reversals, StepSizes.Length - 1) 
                                : 0;
                
                double currentStepSize = (StepSizes != null && StepSizes.Length > 0) 
                                         ? StepSizes[safeIndex] 
                                         : 1.0; 

                if (isCorrect) {
                    consecutiveIncorrect = 0;
                    consecutiveCorrect++;

                    if (consecutiveCorrect >= CorrectRule) {
                        currentIntensity += currentStepSize; // Use dynamic step size
                        consecutiveCorrect = 0;       // Reset counter
                        changed = true;
                        currentDirection = Direction.Up;
                    }
                } else {
                    consecutiveCorrect = 0;
                    consecutiveIncorrect++;

                    if (consecutiveIncorrect >= IncorrectRule) {
                        currentIntensity -= currentStepSize; // Use dynamic step size
                        consecutiveIncorrect = 0;     // Reset counter
                        changed = true;
                        currentDirection = Direction.Down;
                    }
                }

                // Check for a reversal (if the direction changed and it's not the first move)
                if (changed)
                {
                    if (lastDirection != Direction.None && lastDirection != currentDirection)
                    {
                        reversals++;
                        
                        // NEW LOGIC: Record the peak/valley and calculate the moving average
                        reversalPoints.Add(preStepIntensity);
                        
                        // Keep only the last N reversals in the list
                        // CHANGED: 'if' to 'while' so the window can be instantly shrunk on the fly
                        while (reversalPoints.Count > MovingAverageWindow)
                        {
                            reversalPoints.RemoveAt(0);
                        }
                    }
                    lastDirection = currentDirection;
                }
                
                // Clamp the intensity to stay within our defined min/max bounds
                currentIntensity = Math.Max(MinIntensity, Math.Min(MaxIntensity, currentIntensity));
                
                // Calculate current average (returns 0 if no reversals have happened yet)
                double currentAverage = reversalPoints.Count > 0 ? reversalPoints.Average() : 0.0;
                
                return new StaircaseState { 
                    Intensity = currentIntensity, 
                    Reversals = reversals,
                    ReversalAverage = currentAverage // Include average in the output
                };
            })
            // Output initial state before any keys are pressed
            .StartWith(new StaircaseState { 
                Intensity = currentIntensity, 
                Reversals = 0, 
                ReversalAverage = 0.0 
            });
        });
    }
}