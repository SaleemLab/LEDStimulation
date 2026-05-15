using System;
using System.Linq;
using System.Reactive.Linq;
using System.ComponentModel;
using System.Collections.Generic;
using Bonsai;

[Description("Method of Constant Stimuli (MOCS) generator with automatic block randomization.")]
public class MocsProcedure : Transform<object, double>
{
    public MocsProcedure()
    {
        // Default values that appear in the properties panel
        Intensities = new double[] { 2.0, 4.0, 6.0, 8.0, 10.0 };
        Repetitions = 10;
    }

    [Description("The array of stimulus intensities to test.")]
    public double[] Intensities { get; set; }

    [Description("How many times each intensity is placed in the bucket before reshuffling. (Note: The node runs infinitely; this just controls the block size).")]
    public int Repetitions { get; set; }

    public override IObservable<double> Process(IObservable<object> source)
    {
        return Observable.Defer(() => {
            Random rng = new Random();
            List<double> trialQueue = new List<double>();

            // Helper function to refill, shuffle, and pop the next trial
            Func<double> getNext = () => {
                if (trialQueue.Count == 0)
                {
                    // Safety check: Prevent crash if user types 0
                    int safeReps = Repetitions <= 0 ? 1 : Repetitions;

                    // Refill the queue based on properties
                    foreach (var intensity in Intensities)
                    {
                        for (int i = 0; i < safeReps; i++)
                        {
                            trialQueue.Add(intensity);
                        }
                    }
                    
                    // Fisher-Yates Shuffle for true randomization
                    int n = trialQueue.Count;
                    while (n > 1) {
                        n--;
                        int k = rng.Next(n + 1);
                        double value = trialQueue[k];
                        trialQueue[k] = trialQueue[n];
                        trialQueue[n] = value;
                    }
                }

                // Pop the last value
                double nextIntensity = trialQueue[trialQueue.Count - 1];
                trialQueue.RemoveAt(trialQueue.Count - 1);
                return nextIntensity;
            };

            // Generate the very first value at t=0
            double initialIntensity = getNext();

            // Return the stream, starting with the initial value
            return source.Select(_ => getNext())
                         .StartWith(initialIntensity);
        });
    }
}