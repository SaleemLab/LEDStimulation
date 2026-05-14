using Bonsai;
using System;
using System.ComponentModel;
using System.Collections.Generic;
using System.Linq;
using System.Reactive.Linq;
using System.Reactive;

[Combinator]
[Description("")]
[WorkflowElementCategory(ElementCategory.Transform)]
public class RandSampleFromIndexedDistribution
{
    public List<float[]> ProbabilityDistributions {get; set; }
    public int DistributionIndex {get; set; }

    public IObservable<int> Process(IObservable<Unit> source)
    {
        return source.Select(value =>
        {
            float[] thisDistribution; // access probability distribution for stimuli
            int distIndex=DistributionIndex;
            thisDistribution = ProbabilityDistributions[distIndex];
            int nStim = thisDistribution.Length;
            int iStim = 0;

            //Console.WriteLine("dist index: ");
            //Console.WriteLine(distIndex);
            //Console.WriteLine("sample dist: ");
            //Console.WriteLine("[{0}]", string.Join(", ", thisDistribution));


           //for (int i=0; i<nStim; i++) // loop through speed pair
           // {
            //    Console.WriteLine("thisDist" + i + ": " +thisDistribution[i]);
            //}

            var rng = new Random((int)DateTime.Now.Ticks); // random number generator

            int randomStim = rng.Next(1,100); // random number between 1 and 100
            //Console.WriteLine("rng: " + randomStim);
            for (int i=0; i<nStim; i++) // loop through speed pair
            {
                float sum = thisDistribution.Take(i+1).Sum(); // 
                if (randomStim <= sum) // if the sum is more than the random number, choose it
                {
                    iStim = i; // index of the stimulus to show
                    break;
                }
            }

        
            


            return iStim; 
        }
        );
    }
}
