using Bonsai;
using System;
using System.ComponentModel;
using System.Collections.Generic;
using System.Linq;
using System.Reactive.Linq;

[Combinator]
[Description("")]
[WorkflowElementCategory(ElementCategory.Transform)]
public class IncrementStateDependentTrialCounter
{
    public List<int[]> TrialCounterArray {get; set;}
    
    public IObservable<List<int[]>> Process(IObservable<Tuple<int, int>> source)
    {
        

        return source.Select(value => 
        {
            int istate = value.Item1;
            int istim = value.Item2;

            TrialCounterArray[istate][istim]++;

           // Console.WriteLine("state: " + istate);
            //Console.WriteLine("stim: " + istim);
            //Console.WriteLine(TrialCounterArray[istate][istim]);

            return TrialCounterArray;


        });
    }
}
