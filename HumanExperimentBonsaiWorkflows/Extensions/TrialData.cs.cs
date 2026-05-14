using System;
using System.Reactive.Linq;
using System.ComponentModel;
using Bonsai;

// 1. The Master Data Container
public class TrialPayload
{
    // Constructor safely sets defaults
    public TrialPayload()
    {
        Context = "Standing";
        Method = "Staircase";
        StimulusType = "s"; 
        TargetSide = -1;     
        Intensity = 0.0;
        FreqA = 200.0;
        FreqB = 200.0;
        Duration = 1000.0;
        EnvFrequency = 2.0;
        PhaseA = 0.0;
        PhaseB = 0.0;
        ContrastA = 0.5;
        ContrastB = 0.5;
    }

    // Metadata
    public string Context { get; set; }
    public string Method { get; set; }
    public string StimulusType { get; set; }
    
    // Trial Specifics
    public int TargetSide { get; set; }  
    public double Intensity { get; set; } 

    // Hardware Parameters 
    public double FreqA { get; set; } 
    public double FreqB { get; set; } 
    
    // Constants
    public double Duration { get; set; }
    public double EnvFrequency { get; set; } 
    public double PhaseA { get; set; }       
    public double PhaseB { get; set; } 
    public double ContrastA { get; set; }
    public double ContrastB { get; set; }      
}

// 2. The Result Container
public class TrialResult
{
    public TrialPayload OriginalPayload { get; set; }
    public int DiscriminationResponse { get; set; } 
    public int DetectionResponse { get; set; }      
    public bool IsCorrect { get; set; }
}

// 3. The Assembler Node
[Description("Calculates FreqA and FreqB based on Target Side (-1=Left, 1=Right) and the provided Stimulus Intensity (e.g., from Staircase or MOCS).")]
public class CreateTrialPayload : Transform<Tuple<double, int>, TrialPayload> 
{
    public CreateTrialPayload()
    {
        Context = "Standing";
        Method = "Staircase";
        StimulusType = "s";
        BaseFrequency = 200.0;
        Duration = 1000.0;
        EnvFrequency = 2.0;
        ContrastA = 0.5;
        ContrastB = 0.5;
    }

    public string Context { get; set; }
    public string Method { get; set; }
    public string StimulusType { get; set; }
    public double BaseFrequency { get; set; }
    public double Duration { get; set; }
    public double EnvFrequency { get; set; }
    public double ContrastA { get; set; }
    public double ContrastB { get; set; }

    public override IObservable<TrialPayload> Process(IObservable<Tuple<double, int>> source)
    {
        return source.Select(input => {
            double flickerFrequency = input.Item1; 
            int targetSide = input.Item2; 
            
            // Default both sides to the 200Hz baseline
            double freqA = BaseFrequency;
            double freqB = BaseFrequency;

            // -1 is Left, 1 is Right
            if (targetSide == -1) {
                freqA = flickerFrequency;
            } else if (targetSide == 1) {
                freqB = flickerFrequency;
            }

            return new TrialPayload
            {
                Context = this.Context,
                Method = this.Method,
                StimulusType = this.StimulusType,
                TargetSide = targetSide,
                Intensity = flickerFrequency,
                FreqA = freqA,
                FreqB = freqB,
                Duration = this.Duration,
                EnvFrequency = this.EnvFrequency,
                PhaseA = 0, 
                PhaseB = 0,
                ContrastA = this.ContrastA,
                ContrastB = this.ContrastB
            };
        });
    }
}

// 4. The Result Packager Node (MUST BE ITS OWN CLASS)
[Description("Combines the TrialPayload, Discrimination Response, and Detection Response into a single TrialResult object.")]
public class CreateTrialResult : Transform<Tuple<TrialPayload, int, int>, TrialResult>
{
    public override IObservable<TrialResult> Process(IObservable<Tuple<TrialPayload, int, int>> source)
    {
        return source.Select(input => {
            
            TrialPayload payload = input.Item1;
            int discrimResponse = input.Item2;
            int detectResponse = input.Item3;

            // Automatically calculate if they got the Left/Right task correct
            bool isCorrect = (payload.TargetSide == discrimResponse);

            return new TrialResult
            {
                OriginalPayload = payload,
                DiscriminationResponse = discrimResponse,
                DetectionResponse = detectResponse,
                IsCorrect = isCorrect
            };
        });
    }
}