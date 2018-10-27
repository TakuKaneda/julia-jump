"""
Implement perfect foresight with JuMP in parallel
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
@everywhere include("parallel/parallel-perfectforesight.jl")  # load the PerfectForesight function

## SharedArray which are shared by all processors
SampleID = convert(SharedArray,Array(1:NSamples));  # ID of samples
StageCost_vec = convert(SharedArray,Array{Float64}(NSamples,H));  # store cost at each stage
IterationTime = convert(SharedArray,zeros(Float64,(NSamples)));  # store timings

## Generate samples scenarios
# sample_path = SamplePath(TransProb, NSamples);
sample_path = ReadSamplePath("data/test_"* problem_size * "_samples.txt") # if you want to implement with the sample paths

## Prallel for loop
@printf("==== Start Perfect Foresight ====\n")
@parallel for i = 1:NSamples
    RealPath = sample_path[:,:,i]
    solution = Solutions()
    tic()
    StageCost_vec[SampleID[i],:] = PerfectForesight(RealPath, solution);
    IterationTime[SampleID[i]] = toc()
        # @printf(" cost of stage %d sample No.%d:   %5.2f \$\n",t,i,SolutionsArray[i].StageCost[t])
    # end
    @printf("\n====Total cost of sample No.%d:   %5.2f \$====\n\n", SampleID[i],sum(StageCost_vec[SampleID[i],:]))
end
