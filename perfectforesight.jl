### Implement Perfect Foresight policy with JuMP

using DataFrames, JuMP, Gurobi, CSV
include("src/source.jl")

# number of samples
NSamples = 2;

# define solver
solver = GurobiSolver(LogToConsole=0, LogFile="log/train-CeMPC.log")

## Read CSV data
lines_df = CSV.read("data/twolayer-lines.csv")
nodes_df = Read_nodes_csv("data/twolayer-nodes.csv")  # see src/source.jl
generators_df = CSV.read("data/twolayer-generators.csv")

## Read JSON: Network data
LayerNodes, LayerLines = ConvertLayerData2Array("data/two_LayerNodes.json","data/two_LayerLines.json")
HeadNodes, LeafNodes = ConvertHeadLeafNodes2Array("data/two_HeadNodes.json", "data/two_LeafNodes.json")
LeafChildren = ConvertLeafChildren2Array("data/two_LeafChildren.json", LeafNodes,size(nodes_df,1))

## Read JSON: Stochastic Params
PNetDemand = ConvertPNetDemand2Array("data/two_ND.json")
TransProb = ConvertTransProb2Array("data/two_TP.json")
PGenerationMax = ConvertPGenerationCapacity2Array("data/two_PMax.json")
PGenerationMin = ConvertPGenerationCapacity2Array("data/two_PMin.json")
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

## Store Solutions
struct Solutions
    "struct that stores solutions over a sample"
    pflow::Array{Float64,2}
    pgeneration::Array{Float64,2}
    storage::Array{Float64,2}
    batterycharge::Array{Float64,2}
    batterydischarge::Array{Float64,2}
    loadshedding::Array{Float64,2}
    productionshedding::Array{Float64,2}
    p_in::Array{Float64,1}
    p_out::Array{Float64,1}
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
        zeros(Float64,H),
        zeros(Float64,H),
        zeros(Float64,H)
    )
end

function PerfectForesight(RealPath, solutions)
    "
    implement perfect foresight policy
    Args:
        - RealPath: Real path of the 'current' sample
        - solutions: an instance of the struct Solutions
    note: results will be stored in `solutions` thus return is Void
    "

    ## Build model
    m = Model(solver=solver)

    ## Variables
    @variable(m, pflow[1:NLines,1:H])
    @variable(m, pgeneration[1:NGenerators,1:H])
    @variable(m, storage[1:NNodes,1:H] >= 0)
    @variable(m, batterycharge[1:NNodes,1:H] >= 0)
    @variable(m, batterydischarge[1:NNodes,1:H] >= 0)
    @variable(m, loadshedding[1:NNodes,1:H] >= 0)
    @variable(m, productionshedding[1:NNodes,1:H] >= 0)
    @variable(m, p_in[l=1:NLayers,HeadNodes[l],1:H])
    @variable(m, p_out[l=1:NLayers,n in LeafNodes[l], LeafChildren[l,n],1:H])

    ## Objective - minimize cost of generation and load shedding
    @objective(m, Min,
        (sum(MargCost[g]*pgeneration[g,t] for g=1:NGenerators, t=1:H)
        + VOLL * sum(loadshedding))
    );

    ## Constraints
    # dynamics
    @constraint(m, BatteryDynamics_stage1[n = 1:NNodes],
         (storage[n,1] - ini_storage[n]
         - BatteryChargeEfficiency[n] * batterycharge[n,1]
         + batterydischarge[n,1]/BatteryDischargeEfficiency[n]
         == 0)
    );
    @constraint(m, BatteryDynamics[n = 1:NNodes, t=2:H],
        (storage[n,t] - storage[n,t-1]
         - BatteryChargeEfficiency[n] * batterycharge[n,t]
         + batterydischarge[n,t]/BatteryDischargeEfficiency[n]
        == 0)
    );

    # Flow Limits
    @constraint(m, FlowMax[i = 1:NLines, t = 1:H],
        (pflow[i,t] <= SLimit[i])
    );
    @constraint(m, FlowMin[i = 1:NLines, t = 1:H],
        ( - pflow[i,t] <= SLimit[i])
    );

    # Storage Capacity
    @constraint(m, StorageMax[n=1:NNodes, t = 1:H],
        (storage[n,t] <= BatteryCapacity[n])
    );

    # Charging Capacity
    @constraint(m, BatteryChargeMax[n=1:NNodes, t = 1:H],
        (batterycharge[n,t] <= BatteryChargeRate[n])
    );

    # Discharging Capacity
    @constraint(m, BatteryDischargeMax[n=1:NNodes, t = 1:H],
        (batterydischarge[n,t] <= BatteryChargeRate[n])
    );

    # p_in & pflow equality
    @constraint(m, Pin_Flow_equality[l = 1:NLayers, n in HeadNodes[l], t = 1:H],
        (p_in[l,n,t] - pflow[n-1,t] == 0)
    );
    # p_in & p_out equality
    @constraint(m, Pin_Pout_equality[l = 1:NLayers, n in HeadNodes[l], t = 1:H],
        (p_in[l,n,t] - p_out[Node2Layer[Ancestor[n]+1],Ancestor[n]+1,n,t] == 0)
    );

    # Balancing - root node
    @constraint(m, Balance1_rootnode[t = 1:H],
        (sum(pgeneration[g,t] for g = 1:NGenerators)
        + batterydischarge[1,t]+ loadshedding[1,t]
        - productionshedding[1,t] - batterycharge[1,t]
        + sum(pflow[m,t] for m in Children[1])
        == PNetDemand[1,t][RealPath[Node2Layer[1],t]]
        )
    );
    # Balancing - usual nodes
    @constraint(m, Balance[l = 1:NLayers, n in LayerNodes[l], t = 1:H; ~(n==1 || n in HeadNodes[l] || n in LeafNodes[l])],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        - pflow[n-1,t]
        + sum(pflow[m,t] for m in Children[n])
        == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
        )
    );
    # Balancing - head node
    @constraint(m, Balance_headnode[l = 1:NLayers, n in HeadNodes[l], t = 1:H; ~(n in LeafNodes[l])],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        + pflow[n-1,t]
        + sum(pflow[m,t] for m in Children[n])
        == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
        )
    );
    # Balancing - leaf node
    @constraint(m, Balance_leafnode[l=1:NLayers, n in LeafNodes[l], t = 1:H; ~(n in HeadNodes[l])],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        - pflow[n-1,t]
        + sum(pflow[j,t] for j in Children[n] if ~(j+1 in LeafChildren[l,n]))
        - sum(p_out[l,n,j,t] for j in LeafChildren[l,n])
        == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
        )
    );
    # Balancing - head-leaf node
    @constraint(m, Balance_headleafnode[l=1:NLayers, n in LeafNodes[l], t = 1:H; n in HeadNodes[l]],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        + pflow[n-1,t]
        + sum(pflow[j,t] for j in Children[n] if ~(j+1 in LeafChildren[l,n]))
        - sum(p_out[l,n,j,t] for j in LeafChildren[l,n])
        == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
        )
    );

    # Generation Limits
    @constraint(m, GenerationMax[g = 1:NGenerators, t=1:H],
        (pgeneration[g,t] <= PGenerationMax[g,t][RealPath[1,t]])
    );
    @constraint(m, GenerationMin[g = 1:NGenerators, t=1:H],
        ( - pgeneration[g,t] <= - PGenerationMin[g,t][RealPath[1,t]])
    );


    ## Solve
    status = solve(m);
    # println("Objective value: ", getobjectivevalue(m))

    ## Store Results
    solutions.pflow[:,:] = getvalue(pflow);
    solutions.pgeneration[:,:] = getvalue(pgeneration);
    solutions.storage[:,:] = getvalue(storage);
    solutions.batterycharge[:,:] = getvalue(batterycharge);
    solutions.batterydischarge[:,:] = getvalue(batterydischarge);
    solutions.loadshedding[:,:] = getvalue(loadshedding);
    solutions.productionshedding[:,:] = getvalue(productionshedding);
    # solutions.p_in[:] = getvalue(p_in);
    # solutions.p_out[:] = getvalue(p_out);
    solutions.StageCost[:] = [(sum(MargCost[i]*solutions.pgeneration[i,t] for i = 1:NGenerators)
            +VOLL*sum(solutions.loadshedding[:,t])) for t=1:H]
    # solutions.StageCost[:] = [(sum(MargCost[i]*solutions.pgeneration[i,t] for i = 1:NGenerators)+ VOLL * sum(solutions.loadshedding[n,t], for n = 1:NNodes)) for t = 1:H]

    return
end

## Generate samples scenarios
sample_path = SamplePath(TransProb,NSamples);

## Implementation
SolutionsArray = [Solutions() for i=1:NSamples] # array contains Solutions structs
@printf("==== Start Perfect Foresight ====\n")
for i = 1:NSamples
    RealPath = sample_path[:,:,i]
    # for t = 1:H
    PerfectForesight(RealPath, SolutionsArray[i]);
        # @printf(" cost of stage %d sample No.%d:   %5.2f \$\n",t,i,SolutionsArray[i].StageCost[t])
    # end
    @printf("\n====Total cost of sample No.%d:   %5.2f \$====\n\n",i ,sum(SolutionsArray[i].StageCost[:]))
end
