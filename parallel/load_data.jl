"""
Load data from the data folder
"""
# using DataFrames, JuMP, Gurobi, CSV, JSON

## Read CSV data
lines_df = CSV.read("data/" * problem_size * "layer-lines.csv")
nodes_df = Read_nodes_csv("data/" * problem_size * "layer-nodes.csv")  # see src/source.jl
generators_df = CSV.read("data/" * problem_size * "layer-generators.csv")

## Read JSON
PNetDemand = ConvertPNetDemand2Array("data/" * problem_size * "_ND.json")
TransProb = ConvertTransProb2Array("data/" * problem_size * "_TP.json")
PGenerationMax = ConvertPGenerationCapacity2Array("data/" * problem_size * "_PMax.json")
PGenerationMin = ConvertPGenerationCapacity2Array("data/" * problem_size * "_PMin.json")

## Problem Parameters
# generators
Generators = generators_df[:GeneratorID]
MargCost = generators_df[:MargCost]
# lines
Lines = lines_df[:LineID]
SLimit = lines_df[:SLimit]
# nodes
Nodes = nodes_df[:NodeID]
BatteryCapacity = nodes_df[:BatteryCapacity]
BatteryChargeEfficiency = nodes_df[:BatteryChargeEfficiency]
BatteryDischargeEfficiency = nodes_df[:BatteryDischargeEfficiency]
BatteryChargeRate = nodes_df[:BatteryChargeRate]
ini_storage = nodes_df[:ini_storage]
Ancestor = nodes_df[:Ancestor]
Children = nodes_df[:Children]
Node2Layer = nodes_df[:Node2Layer]
# other parameters
VOLL = 5000;
H = size(PNetDemand,2);
T = 1:H;
NLayers = size(TransProb,1)
NLines = nrow(lines_df)
NNodes = nrow(nodes_df)
NGenerators = nrow(generators_df)
NLattice = [size(PNetDemand[1,t],1) for t = 1:H]
# Policy specific parameters: expected value of stochastic parameters
PNetDemand_fix = Array{Float64}(NNodes, H);
PGenerationMax_fix = Array{Float64}(NGenerators, H);
PGenerationMin_fix = Array{Float64}(NGenerators, H);

## Store Solutions
immutable Solutions # change to immutable so as to be converted to SharedArray
    "struct that stores solutions over a sample"
    pflow::Array{Float64,2}
    pgeneration::Array{Float64,2}
    storage::Array{Float64,2}
    batterycharge::Array{Float64,2}
    batterydischarge::Array{Float64,2}
    loadshedding::Array{Float64,2}
    productionshedding::Array{Float64,2}
    # p_in::Array{Float64,1}
    # p_out::Array{Float64,1}
    StageCost::Array{Float64,1}

    # constructor
    # maybe there is a better way to assign default values
    Solutions() = new(
        zeros(Float64,(NLines, H)),
        zeros(Float64,(NGenerators, H)),
        zeros(Float64,(NNodes, H)),
        zeros(Float64,(NNodes, H)),
        zeros(Float64,(NNodes, H)),
        zeros(Float64,(NNodes, H)),
        zeros(Float64,(NNodes, H)),
        # zeros(Float64,H),
        # zeros(Float64,H),
        zeros(Float64,H)
    )
end
