### Implement Perfect Foresight policy with JuMP
using DataFrames, JuMP, Gurobi, CSV
include("src/source.jl")

## choose the problem size
# problem_size = {"two", "multi"}
problem_size = "multi"

## number of samples
NSamples = 2;

## define solver
solver = GurobiSolver(LogToConsole=0, LogFile="log/train-PerfectForesight.log")

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
    # @variable(m, p_in[1:H])
    # @variable(m, p_out[1:H])

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
    # @constraint(m, Pin_Flow_equality[t = 1:H],
    #     (p_in[t] - pflow[8,t] == 0)
    # );
    # # p_in & p_out equality
    # @constraint(m, Pin_Pout_equality[t = 1:H],
    #     (p_in[t] - p_out[t] == 0)
    # );

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
    @constraint(m, Balance[n = 1:NNodes, t = 1:H ; n!=0+1 #=&& n!=3+1 && n!=8+1=#],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        - pflow[n-1,t]
        + sum(pflow[m,t] for m in Children[n])
        == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
        )
    );
    # Balancing - head node
    # @constraint(m, Balance_headnode[n in [8+1], t = 1:H],
    #     (batterydischarge[n,t]+ loadshedding[n,t]
    #     - productionshedding[n,t]- batterycharge[n,t]
    #     + pflow[n-1,t]
    #     + sum(pflow[m,t] for m in Children[n])
    #     == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
    #     )
    # );
    # # Balancing - leaf node
    # @constraint(m, Balance_leafnode[n in [3+1], t = 1:H],
    #     (batterydischarge[n,t]+ loadshedding[n,t]
    #     - productionshedding[n,t]- batterycharge[n,t]
    #     - pflow[n-1,t]
    #     + sum(pflow[m,t] for m in [4])
    #     - p_out[t]
    #     == PNetDemand[n,t][RealPath[Node2Layer[n],t]]
    #     )
    # );

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
    tic()
    PerfectForesight(RealPath, SolutionsArray[i]);
    toc()
        # @printf(" cost of stage %d sample No.%d:   %5.2f \$\n",t,i,SolutionsArray[i].StageCost[t])
    # end
    @printf("\n====Total cost of sample No.%d:   %5.2f \$====\n\n",i ,sum(SolutionsArray[i].StageCost[:]))
end
