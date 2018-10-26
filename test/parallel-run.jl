NSamples = 10;
NProcessors = 2;

if nworkers() == 1
    addprocs(NProcessors);
end
## load the problem
@everywhere include("parallel.jl")
SampleID = convert(SharedArray,Array(1:NSamples));
sol_vec = convert(SharedArray,Array{Float64}(NSamples));
## problem parameters
h_vec = rand(1:30,NSamples,3)
##
@parallel for i = 1:NSamples
    println("Current Sample ID: ",SampleID[i])
    println("Start Solve ID: ",SampleID[i])
    sol_vec[SampleID[i]] = sample_problem(h_vec[SampleID[i],:])
    println("End Solve ID: ",SampleID[i])
end
