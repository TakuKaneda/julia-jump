"""
Implement ce-mpc with JuMP in parallel
"""
## number of samples
NSamples = 3;

## number of processors
NProcessors = 2;
addprocs(NProcessors)  # add processors

## load global modules, problem data and MPC function
@everywhere problem_size = "two" # two or multi: problem name
@everywhere using DataFrames, JuMP, Gurobi, CSV, JSON
@everywhere include("src/source.jl")  # functions used to load problem data
@everywhere include("parallel/load_data.jl")  # load problem data
@everywhere include("parallel/parallel-cempc.jl")  # load the MPC function

## SharedArray which are shared by all processors
SampleID = convert(SharedArray,Array(1:NSamples));  # ID of samples
StageCost_vec = convert(SharedArray,Array{Float64}(NSamples,H));  # store cost at each stage
IterationTime = convert(SharedArray,zeros(Float64,(NSamples,H)));  # store timings

## Generate samples scenarios
# sample_path = SamplePath(TransProb, NSamples);
sample_path = ReadSamplePath("data/test_"* problem_size * "_samples.txt") # if you want to implement with the sample paths

## Prallel for loop
@printf("==== Start certainty-equivalent MPC ====\n")
@parallel for i = 1:NSamples
    RealPath = sample_path[:,:,SampleID[i]]
    solution = Solutions()
    for t = 1:H
        tic();
        StageCost_vec[SampleID[i],t] =  CeMPC(t, RealPath, solution);
        IterationTime[SampleID[i],t] = toc();
        @printf(" cost of stage %d sample No.%d:   %5.2f \$\n",t,SampleID[i],StageCost_vec[SampleID[i],t])
    end
    @printf("\n====Total cost of sample No.%d:   %5.2f \$====\n\n", SampleID[i],sum(StageCost_vec[SampleID[i],:]))
end
