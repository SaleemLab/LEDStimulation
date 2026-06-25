using System;
using System.Linq;
using System.Reactive.Linq;
using System.ComponentModel;
using System.Collections.Generic;
using Bonsai;

[Description("Method of Constant Stimuli (MOCS) generator. Outputs randomized Tuple (Intensity as Item1, Side as Item2).")]
public class MocsProcedure : Transform<object, Tuple<double, int>>
{
    public MocsProcedure()
    {
        // Default values that appear in the properties panel
        Intensities = new double[] { 2.0, 4.0, 6.0, 8.0, 10.0 };
        Repetitions = 10;
    }

    [Description("The array of stimulus intensities to test.")]
    public double[] Intensities { get; set; }

    [Description("How many times each intensity (paired with a Left and Right) is placed in the bucket before reshuffling.")]
    public int Repetitions { get; set; }

    public override IObservable<Tuple<double, int>> Process(IObservable<object> source)
    {
        return Observable.Defer(() => {
            Random rng = new Random();
            
            List<Tuple<double, int>> trialQueue = new List<Tuple<double, int>>();

            // Helper function to refill, shuffle, and pop the next trial
            Func<Tuple<double, int>> getNext = () => {
                if (trialQueue.Count == 0)
                {
                    // Safety check: Prevent crash if user types 0
                    int safeReps = Repetitions <= 0 ? 1 : Repetitions;

                    // Refill the queue: For each intensity, add exactly 1 and -1 per repetition
                    foreach (var intensity in Intensities)
                    {
                        for (int i = 0; i < safeReps; i++)
                        {
                            trialQueue.Add(Tuple.Create(intensity, -1));
                            trialQueue.Add(Tuple.Create(intensity, 1));
                        }
                    }
                    
                    // Fisher-Yates Shuffle for true randomization
                    int n = trialQueue.Count;
                    while (n > 1) {
                        n--;
                        int k = rng.Next(n + 1);
                        Tuple<double, int> value = trialQueue[k];
                        trialQueue[k] = trialQueue[n];
                        trialQueue[n] = value;
                    }
                }

                // Pop the last value
                Tuple<double, int> nextTrial = trialQueue[trialQueue.Count - 1];
                trialQueue.RemoveAt(trialQueue.Count - 1);
                return nextTrial;
            };

            // Generate the very first value at t=0
            Tuple<double, int> initialTrial = getNext();

            // Return the stream, starting with the initial value
            return source.Select(_ => getNext())
                         .StartWith(initialTrial);
        });
    }
}