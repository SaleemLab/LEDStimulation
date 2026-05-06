using System;
using System.Reactive.Linq;
using System.ComponentModel;
using Bonsai;

// We define a custom object to hold both the Intensity and the Reversal count
public class StaircaseState
{
    public double Intensity { get; set; }
    public int Reversals { get; set; }
}

[Description("Adaptive Staircase tracking intensity and reversals.")]
public class StaircaseProcedure : Transform<bool, StaircaseState>
{
    public StaircaseProcedure()
    {
        InitialIntensity = 10.0;
        StepSize = 1.0;
        MinIntensity = 0.0;
        MaxIntensity = 100.0;
        CorrectRule = 2;   // Formerly DownRule
        IncorrectRule = 1; // Formerly UpRule
    }

    [Description("The initial stimulus intensity before any responses are made.")]
    public double InitialIntensity { get; set; }

    [Description("The amount to change the intensity.")]
    public double StepSize { get; set; }

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

            // source is our incoming stream of booleans (true = correct, false = incorrect)
            return source.Select(isCorrect => {
                bool changed = false;
                Direction currentDirection = Direction.None;

                if (isCorrect) {
                    consecutiveIncorrect = 0;
                    consecutiveCorrect++;

                    if (consecutiveCorrect >= CorrectRule) {
                        currentIntensity += StepSize; // Make value higher
                        consecutiveCorrect = 0;       // Reset counter
                        changed = true;
                        currentDirection = Direction.Up;
                    }
                } else {
                    consecutiveCorrect = 0;
                    consecutiveIncorrect++;

                    if (consecutiveIncorrect >= IncorrectRule) {
                        currentIntensity -= StepSize; // Make value lower
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
                    }
                    lastDirection = currentDirection;
                }
                
                // Clamp the intensity to stay within our defined min/max bounds
                currentIntensity = Math.Max(MinIntensity, Math.Min(MaxIntensity, currentIntensity));
                
                return new StaircaseState { 
                    Intensity = currentIntensity, 
                    Reversals = reversals 
                };
            })
            // Output initial state before any keys are pressed
            .StartWith(new StaircaseState { Intensity = currentIntensity, Reversals = 0 });
        });
    }
}