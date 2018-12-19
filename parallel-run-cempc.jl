# Implement ce-mpc with JuMP in parallel

## number of samples
NSamples = 4;

## number of processors
NProcessors = 4;
addprocs(NProcessors)  # add processors

## load global modules, problem data and MPC function
@everywhere problem_size = "multi" # two or multi: problem name
@everywhere using DataFrames, JuMP, Gurobi, CSV, JSON
@everywhere include("src/source.jl")  # functions used to load problem data
@everywhere include("parallel/load_data.jl")  # load problem data
@everywhere include("parallel/parallel-cempc.jl")  # load the MPC function

## Prallel implementation by pmap
@printf("==== Start certainty-equivalent MPC ====\n")
tic()
SolutionsArray =  pmap(CeMPC_all, 1:NSamples)
toc()
