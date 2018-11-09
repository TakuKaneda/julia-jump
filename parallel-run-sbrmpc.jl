"""
Implement sbr-mpc with JuMP in parallel
"""
## number of samples
NSamples = 4;

## number of processors
NProcessors = 2;
addprocs(NProcessors)  # add processors

########################
# Hyperparameters of the algorithm: SET AS YOU WANT
@everywhere NScenarios = 5; # number of generated scenarios at each stage
@everywhere DiscountFactor = 0.9;
########################

## load global modules, problem data and MPC function
@everywhere problem_size = "multi" # two or multi: problem name
@everywhere using DataFrames, JuMP, Gurobi, CSV, JSON
@everywhere include("src/source.jl")  # functions used to load problem data
@everywhere include("parallel/load_data.jl")  # load problem data
@everywhere include("parallel/parallel-sbrmpc.jl")  # load the MPC function

## Prallel implementation by pmap
@printf("==== Start scenario-based robust MPC ====\n")
tic()
SolutionsArray =  pmap(SbrMPC_all, 1:NSamples)
toc()
