### Implement sbr-mpc with JuMP
using DataFrames, JuMP, Gurobi, CSV, JSON
include("src/source.jl")

## choose the problem size
# problem_size = {"two", "multi"}
problem_size = "multi"

########################
# Hyperparameters of the algorithm: SET AS YOU WANT
NScenarios = 5;
DiscountFactor = 0.9;
########################

# number of samples
NSamples = 1;

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


function ScenarioGenerator(TimeChoice, RealPath)
    "
    Generate scenarios for a given stage and outcome via Unform Sampling
    Args:
        - TimeChoice: 'current' stage
        - RealPath[l, t]: Real path of the 'current' sample
    Output:
        - SbrScenarios[l,t,s]: outcome of layer l at stage t for s-th scenario
    "
    SbrScenarios = Array{Int64}(NLayers, H, NScenarios)
    for s=1:NScenarios
        SbrScenarios[:,1:TimeChoice,s] = RealPath[:,1:TimeChoice] # assignment
    end

    if TimeChoice == H
        return SbrScenarios
    else
        for u= TimeChoice+1:H
            SbrScenarios[:,u,:] = rand(1:NLattice[u],NLayers,1,NScenarios)
        end
        return SbrScenarios
    end
end


function SbrMPC(TimeChoice, RealPath, solutions)
    "
    implement sbr-mpc (scenario-based robust MPC)
    Args:
        - TimeChoice: 'current' stage
        - RealPath: Real path of the 'current' sample
        - solutions: an instance of the struct Solutions
    note: results will be stored in `solutions` thus return is Void
    "
    ## Generate scenarios
    SbrScenarios = ScenarioGenerator(TimeChoice, RealPath);

    ## Build model
    m = Model(solver=solver)

    ## Variables
    @variable(m, pflow[1:NLines, TimeChoice:H, 1:NScenarios])
    @variable(m, pgeneration[1:NGenerators, TimeChoice:H, 1:NScenarios])
    @variable(m, storage[1:NNodes, TimeChoice:H, 1:NScenarios] >= 0)
    @variable(m, batterycharge[1:NNodes, TimeChoice:H, 1:NScenarios] >= 0)
    @variable(m, batterydischarge[1:NNodes, TimeChoice:H, 1:NScenarios] >= 0)
    @variable(m, loadshedding[1:NNodes, TimeChoice:H, 1:NScenarios] >= 0)
    @variable(m, productionshedding[1:NNodes, TimeChoice:H, 1:NScenarios] >= 0)
    # @variable(m, p_in[TimeChoice:H, 1:NScenarios])
    # @variable(m, p_out[TimeChoice:H, 1:NScenarios])
    @variable(m, var_cost)

    ## Objective - minimize cost of generation and load shedding till the Horizon
    @objective(m, Min,var_cost);

    ## Constraints
    # upper bound on the cost (which we want to minimize)
    @constraint(m, CostBound[s=1:NScenarios],
        (sum(DiscountFactor^(u-TimeChoice)*
                (sum(MargCost[i]*pgeneration[i,u,s] for i in 1:NGenerators)
                + VOLL * sum(loadshedding[:,u,s])) for u in TimeChoice:H
        ) <= var_cost)
    );

    # dynamics
    if TimeChoice == 1
        # if stage 1
        @constraint(m, BatteryDynamics_stage1[n=1:NNodes,s=1:NScenarios],
             (storage[n,1,s] - ini_storage[n]
             - BatteryChargeEfficiency[n] * batterycharge[n,1,s]
             + batterydischarge[n,1,s]/BatteryDischargeEfficiency[n]
             == 0)
        );
    else
        # current stage
        @constraint(m, BatteryDynamics_current[n=1:NNodes,s=1:NScenarios],
            (storage[n,TimeChoice,s] - solutions.storage[n,TimeChoice-1]  # Unique for a sample
             - BatteryChargeEfficiency[n] * batterycharge[n,TimeChoice,s]
             + batterydischarge[n,TimeChoice,s]/BatteryDischargeEfficiency[n]
            == 0)
        );
    end
    # future stages
    if TimeChoice < H
        @constraint(m, BatteryDynamics_future[n=1:NNodes,u=TimeChoice+1:H,s=1:NScenarios],
            (storage[n,u,s] - storage[n,u-1,s]
             - BatteryChargeEfficiency[n] * batterycharge[n,u,s]
             + batterydischarge[n,u,s]/BatteryDischargeEfficiency[n]
            == 0)
        );
    end

    # Flow Limits
    @constraint(m, FlowMax[i=1:NLines, u = TimeChoice:H, s=1:NScenarios],
        (pflow[i,u,s] <= SLimit[i])
    );
    @constraint(m, FlowMin[i=1:NLines, u = TimeChoice:H, s=1:NScenarios],
        ( - pflow[i,u,s] <= SLimit[i])
    );

    # Storage Capacity
    @constraint(m, StorageMax[n=1:NNodes, u = TimeChoice:H, s=1:NScenarios],
        (storage[n,u,s] <= BatteryCapacity[n])
    );

    # Charging Capacity
    @constraint(m, BatteryChargeMax[n=1:NNodes, u = TimeChoice:H, s=1:NScenarios],
        (batterycharge[n,u,s] <= BatteryChargeRate[n])
    );

    # Discharging Capacity
    @constraint(m, BatteryDischargeMax[n=1:NNodes, u = TimeChoice:H, s=1:NScenarios],
        (batterydischarge[n,u,s] <= BatteryChargeRate[n])
    );

    # # p_in & pflow equality
    # @constraint(m, Pin_Flow_equality[u = TimeChoice:H, s=1:NScenarios],
    #     (p_in[u,s] - pflow[8,u,s] == 0)
    # );
    # # p_in & p_out equality
    # @constraint(m, Pin_Pout_equality[u = TimeChoice:H, s=1:NScenarios],
    #     (p_in[u,s] - p_out[u,s] == 0)
    # );

    # Balancing - root node
    @constraint(m, Balance1_rootnode[u = TimeChoice:H, s=1:NScenarios],
        (sum(pgeneration[g,u,s] for g = 1:NGenerators)
        + batterydischarge[1,u,s]+ loadshedding[1,u,s]
        - productionshedding[1,u,s] - batterycharge[1,u,s]
        + sum(pflow[m,u,s] for m in Children[1])
        == PNetDemand[1,u][SbrScenarios[Node2Layer[1],u,s]])
    );
    # Balancing - usual node
    @constraint(m, Balance[n = 1:NNodes, u = TimeChoice:H, s=1:NScenarios;
    n!=0+1 #=&& n!=3+1 && n!=8+1=#],
        (batterydischarge[n,u,s]+ loadshedding[n,u,s]
        - productionshedding[n,u,s]- batterycharge[n,u,s]
        - pflow[n-1,u,s]
        + sum(pflow[m,u,s] for m in Children[n])
        == PNetDemand[n,u][SbrScenarios[Node2Layer[n],u,s]])
    );
    # # Balancing - head node
    # @constraint(m, Balance_headnode[n in [8+1], u = TimeChoice:H, s=1:NScenarios],
    #     (batterydischarge[n,u,s]+ loadshedding[n,u,s]
    #     - productionshedding[n,u,s]- batterycharge[n,u,s]
    #     + pflow[n-1,u,s]
    #     + sum(pflow[m,u,s] for m in Children[n])
    #     == PNetDemand[n,u][SbrScenarios[Node2Layer[n],u,s]])
    # );
    # # Balancing - leaf node
    # @constraint(m, Balance_leafnode[n in [3+1], u = TimeChoice:H, s=1:NScenarios],
    #     (batterydischarge[n,u,s]+ loadshedding[n,u,s]
    #     - productionshedding[n,u,s]- batterycharge[n,u,s]
    #     - pflow[n-1,u,s]
    #     + sum(pflow[m,u,s] for m in [4])
    #     - p_out[u,s]
    #     == PNetDemand[n,u][SbrScenarios[Node2Layer[n],u,s]])
    # );

    # Generation Limits
    @constraint(m, GenerationMax[g = 1:NGenerators, u = TimeChoice:H, s=1:NScenarios],
        (pgeneration[g,u,s] <= PGenerationMax[g,u][SbrScenarios[1,u,s]])
    );
    @constraint(m, GenerationMin[g = 1:NGenerators, u = TimeChoice:H, s=1:NScenarios],
        ( - pgeneration[g,u,s] <= - PGenerationMin[g,u][SbrScenarios[1,u,s]])
    );

    # Equality of variables at TimeChoice to be executed
    # note: ignore p_in and p_out here
    #       (since they are guaranteed by Pin_Flow_equality and Pin_Pout_equality)
    @constraint(m, Equal_pflow[i = 1:NLines, s=1:NScenarios],
        (pflow[i,TimeChoice,s] == sum(pflow[i,TimeChoice,:]) / NScenarios)
    );
    @constraint(m, Equal_pgeneration[g = 1:NGenerators, s=1:NScenarios],
        (pgeneration[g,TimeChoice,s] == sum(pgeneration[g,TimeChoice,:]) / NScenarios )
    );
    @constraint(m,Equal_storage[n=1:NNodes, s=1:NScenarios],
        (storage[n,TimeChoice,s] == sum(storage[n,TimeChoice,:]) / NScenarios)
    );
    @constraint(m,Equal_batterycharge[n=1:NNodes, s=1:NScenarios],
        (batterycharge[n,TimeChoice,s] == sum(batterycharge[n,TimeChoice,:]) / NScenarios)
    );
    @constraint(m,Equal_batterydischarge[n=1:NNodes, s=1:NScenarios],
        (batterydischarge[n,TimeChoice,s] == sum(batterydischarge[n,TimeChoice,:]) / NScenarios)
    );
    @constraint(m,Equal_loadshedding[n=1:NNodes, s=1:NScenarios],
        (loadshedding[n,TimeChoice,s] == sum(loadshedding[n,TimeChoice,:]) / NScenarios)
    );
    @constraint(m,Equal_productionshedding[n=1:NNodes, s=1:NScenarios],
        (productionshedding[n,TimeChoice,s] == sum(productionshedding[n,TimeChoice,:]) / NScenarios)
    );

    ## Solve
    @time status = solve(m);

    ## Store Results
    solutions.pflow[:,TimeChoice] = getvalue(pflow[:,TimeChoice,1])
    solutions.pgeneration[:,TimeChoice] = getvalue(pgeneration[:,TimeChoice,1])
    solutions.storage[:,TimeChoice] = getvalue(storage[:,TimeChoice,1])
    solutions.batterycharge[:,TimeChoice] = getvalue(batterycharge[:,TimeChoice,1])
    solutions.batterydischarge[:,TimeChoice] = getvalue(batterydischarge[:,TimeChoice,1])
    solutions.loadshedding[:,TimeChoice] = getvalue(loadshedding[:,TimeChoice,1])
    solutions.productionshedding[:,TimeChoice] = getvalue(productionshedding[:,TimeChoice,1])
    # solutions.p_in[TimeChoice] = getvalue(p_in[TimeChoice,1])
    # solutions.p_out[TimeChoice] = getvalue(p_out[TimeChoice,1])
    solutions.StageCost[TimeChoice] = (sum(MargCost[i]*solutions.pgeneration[i,TimeChoice] for i in 1:NGenerators)
                    + VOLL * sum(solutions.loadshedding[:,TimeChoice]))

    return
end

## Generate samples scenarios
sample_path = SamplePath(TransProb,NSamples);

## Implementation
SolutionsArray = [Solutions() for i=1:NSamples] # array contains Solutions structs
@printf("==== Start scenario-based robust MPC ====\n")
for i = 1:NSamples
    RealPath = sample_path[:,:,i]
    for t = 1:H
        tic()
        SbrMPC(t, RealPath, SolutionsArray[i]);
        toc()
        @printf(" cost of stage %d sample No.%d:   %5.2f \$\n",t,i,SolutionsArray[i].StageCost[t])
    end
    @printf("\n====Total cost of sample No.%d:   %5.2f \$====\n\n",i ,sum(SolutionsArray[i].StageCost[:]))
end
