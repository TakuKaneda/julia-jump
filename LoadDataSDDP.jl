using DataFrames, CSV
include("src/source.jl")

## Read CSV data
lines_df = CSV.read("data/twolayer-lines.csv")
nodes_df = Read_nodes_csv("data/twolayer-nodes.csv")  # see src/source.jl
generators_df = CSV.read("data/twolayer-generators.csv")

## Read JSON: Network data
LayerNodes, LayerLines = ConvertLayerData2Array("data/two_LayerNodes.json", "data/two_LayerLines.json")
HeadNodes, LeafNodes = ConvertHeadLeafNodes2Array("data/two_HeadNodes.json", "data/two_LeafNodes.json")
LeafChildren = ConvertLeafChildren2Array("data/two_LeafChildren.json", LeafNodes,size(nodes_df,1))

## Read JSON: Stochastic Params
PNetDemand = ConvertPNetDemand2Array("data/two_ND.json")
TransProb = ConvertTransProb2Array("data/two_TP.json")
PGenerationMax = ConvertPGenerationCapacity2Array("data/two_PMax.json")
PGenerationMin = ConvertPGenerationCapacity2Array("data/two_PMin.json")

## Read txt: Stochastic Params
p_in_data = readdlm("SDDP_data/p_in_lattice.txt")
p_out_data = readdlm("SDDP_data/p_out_lattice.txt")
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
NLayers = size(TransProb,1);
NLines = nrow(lines_df);
NNodes = nrow(nodes_df);
NGenerators = nrow(generators_df);
NLattice = [size(PNetDemand[1,t],1) for t = 1:H];

# Create sets
UsualNodeSet = [];
HeadNodeSet = [];
LeafNodeSet = [];
HeadLeafNodeSet = [];
ChildrenNodesMinusLeafChildren = Dict{Array{Int64,1},Int64}();
for LayerChoice = 1:NLayers
    push!(UsualNodeSet, setdiff(LayerNodes[LayerChoice], union( [1], HeadNodes[LayerChoice], LeafNodes[LayerChoice] )) );
    push!(HeadNodeSet, setdiff(HeadNodes[LayerChoice], LeafNodes[LayerChoice]) );
    push!(LeafNodeSet, setdiff(LeafNodes[LayerChoice], HeadNodes[LayerChoice]) );
    push!(HeadLeafNodeSet, intersect(LeafNodes[LayerChoice], HeadNodes[LayerChoice]) );
    
    for n in LeafNodeSet[LayerChoice]
        for m in Children[n]
            if issubset([m+1], LeafChildren[LayerChoice,n]) == false
                 ChildrenNodesMinusLeafChildren[[LayerChoice,n]] = m;
            end
        end
    end
end
