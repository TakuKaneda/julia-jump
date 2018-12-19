# Implement perfect foresight with JuMP in parallel

## number of samples
NSamples = 10;

## number of processors
NProcessors = 4;
addprocs(NProcessors)  # add processors

## load global modules, problem data and MPC function
@everywhere problem_size = "two" # two or multi: problem name
@everywhere using DataFrames, JuMP, Gurobi, CSV, JSON
@everywhere include("src/source.jl")  # functions used to load problem data
@everywhere include("parallel/load_data.jl")  # load problem data
@everywhere include("parallel/parallel-perfectforesight.jl")  # load the PerfectForesight function

## Prallel implementation by pmap
@printf("==== Start Perfect Foresight ====\n")
tic()
SolutionsArray =  pmap(PerfectForesight_all, 1:NSamples)
toc()
