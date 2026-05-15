using System;
using System.Linq;
using System.Collections.Generic;
using System.Reactive.Linq;
using System.ComponentModel;
using Bonsai;

[Description("Aggregates TrialResults and emits a Tuple containing two arrays: (Frequencies[], pCorrect[]).")]
public class PsychometricTracker : Transform<TrialResult, Tuple<double[], double[]>>
{
    public override IObservable<Tuple<double[], double[]>> Process(IObservable<TrialResult> source)
    {
        return Observable.Defer(() => {
            
            // Internal memory dictionary mapping Intensity -> Tuple(TotalTrials, CorrectTrials)
            var data = new Dictionary<double, Tuple<int, int>>();

            return source.Select(result => {
                
                // Extract what we need from the incoming TrialResult
                double intensity = result.OriginalPayload.Intensity;
                bool isCorrect = result.IsCorrect;

                // Initialize the bin if this is the first time seeing this intensity
                if (!data.ContainsKey(intensity)) {
                    data[intensity] = Tuple.Create(0, 0);
                }

                // Update the counts for this specific intensity bin
                var current = data[intensity];
                data[intensity] = Tuple.Create(current.Item1 + 1, current.Item2 + (isCorrect ? 1 : 0));

                // Sort the dictionary so the curve is ordered left-to-right
                var sortedData = data.OrderBy(kvp => kvp.Key).ToList();

                // Generate the two distinct arrays
                double[] frequencies = sortedData.Select(kvp => kvp.Key).ToArray();
                double[] pCorrects = sortedData.Select(kvp => (double)kvp.Value.Item2 / kvp.Value.Item1).ToArray();

                // Output the Tuple of arrays
                return Tuple.Create(frequencies, pCorrects);
            });
        });
    }
}