### Implement ce-mpc with JuMP
using DataFrames, JuMP, Gurobi, CSV, JSON
include("test-source.jl")
solver = GurobiSolver(Presolve=0, LogToConsole=0, LogFile="log/train-CeMPC.log")

# number of samples
NSamples = 2;

## Read CSV data
lines_df = CSV.read("data/twolayer-lines.csv")
nodes_df = CSV.read("data/twolayer-nodes.csv")
generators_df = CSV.read("data/twolayer-generators.csv")

## Read JSON
PNetDemand = ConvertPNetDemand2Array("data/two_ND.json")
TransProb = ConvertTransProb2Array("data/two_TP.json")
PGenerationMax = ConvertPGenerationCapacity2Array("data/two_PMax.json")
PGenerationMin = ConvertPGenerationCapacity2Array("data/two_PMin.json")

## Data clearning
# convert Children to Array
s = []
for i = 1:nrow(nodes_df)
    if typeof(nodes_df[i,:Children]) == Missing
        push!(s,[])
    else
        push!(s,map(parse,split(nodes_df[i,:Children])))
    end
end
nodes_df[:Children] = s  # assign
## check DataFrames
# print(nodes_df)
# print(lines_df)
# print(generators_df)
# print(PNetDemand_df)

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
NLattice = Array{Int64}(H)
for t=1:H
    NLattice[t] = size(PNetDemand[1,t],1)
end
## Store Solutions
pflow_Solution = Array{Float64}(NLines, H, NSamples)
pgeneration_Solution = Array{Float64}(NGenerators, H, NSamples)
storage_Solution = Array{Float64}(NNodes, H, NSamples)
batterycharge_Solution = Array{Float64}(NNodes, H, NSamples)
batterydischarge_Solution = Array{Float64}(NNodes, H, NSamples)
loadshedding_Solution = Array{Float64}(NNodes, H, NSamples)
productionshedding_Solution = Array{Float64}(NNodes, H, NSamples)
p_in_Solution = Array{Float64}(H, NSamples)
p_out_Solution = Array{Float64}(H, NSamples)


## Compute expected value of stochastic parameters
PNetDemand_fix = Array{Float64}(NNodes, H);
PGenerationMax_fix = Array{Float64}(NGenerators, H);
PGenerationMin_fix = Array{Float64}(NGenerators, H);
##
function ComputeExpectedParameters(TimeChoice,ScenarioChoice)
    "
    Compute expected value of stochastic parameters of the problem, i.e.
    * PNetDemand_fix
    * PGenerationMax_fix
    * PGenerationMin_fix
    Args:
        - TimeChoice: 'current' stage
        - ScenarioChoice[l]: vector of outcomes for layer l
    "
    # assign current observation
    PNetDemand_fix[:,TimeChoice] =  [PNetDemand[n,TimeChoice][ScenarioChoice[Node2Layer[n]]] for n =1:NNodes];
    PGenerationMax_fix[:,TimeChoice] = [PGenerationMax[g,TimeChoice][ScenarioChoice[1]] for g=1:NGenerators];
    PGenerationMin_fix[:,TimeChoice] = [PGenerationMin[g,TimeChoice][ScenarioChoice[1]] for g=1:NGenerators];
    prob = Array{Array}(NLayers,H);
    if TimeChoice == 1
        # compute prob of landing future outcomes
        for l=1:NLayers
            prob[l,1] = [1.0];
            for u=2:H
                prob[l,u] = []
                for k=1:NLattice[u]
                    prob[l,u] = push!(prob[l,u], sum(TransProb[l,u-1][j,k]*prob[l,u-1][j] for j=1:NLattice[u-1]))
                end
                prob[l,u] = convert(Array{Float64,1},prob[l,u])
            end
        end
        # compute ExpValues
        for u=2:H
            PNetDemand_fix[:,u] = [prob[Node2Layer[n],u]'*PNetDemand[n,u] for n=1:NNodes];
            PGenerationMax_fix[:,u] = [prob[1,u]'*PGenerationMax[g,u] for g=1:NGenerators];
            PGenerationMin_fix[:,u] = [prob[1,u]'*PGenerationMin[g,u] for g=1:NGenerators];
        end
    elseif TimeChoice <= H - 1
        for l=1:NLayers
            # prob[l,TimeChoice+1] = TransProb[l,TimeChoice][sample_path[l,TimeChoice,SampleChoice],:]
            prob[l,TimeChoice+1] = TransProb[l,TimeChoice][ScenarioChoice[l],:]
            for u=TimeChoice+2:H
                prob[l,u] = []
                for k=1:NLattice[u]
                    prob[l,u] = push!(prob[l,u], sum(TransProb[l,u-1][j,k]*prob[l,u-1][j] for j=1:NLattice[u-1]))
                end
                prob[l,u] = convert(Array{Float64,1},prob[l,u])
            end
        end
        for u=TimeChoice+1:H
            PNetDemand_fix[:,u] = [prob[Node2Layer[n],u]'*PNetDemand[n,u] for n=1:NNodes];
            PGenerationMax_fix[:,u] = [prob[1,u]'*PGenerationMax[g,u] for g=1:NGenerators];
            PGenerationMin_fix[:,u] = [prob[1,u]'*PGenerationMin[g,u] for g=1:NGenerators];
        end
    end
    return PNetDemand_fix, PGenerationMax_fix, PGenerationMin_fix
end
## implement certainty-equivalent Model Predictive Control
function CeMPC(TimeChoice,ScenarioChoice,SampleChoice)
    ## Compute expected value of stochastic parameters
    PNetDemand_fix, PGenerationMax_fix, PGenerationMin_fix = ComputeExpectedParameters(TimeChoice,ScenarioChoice)

    ## Build model
    m = Model(solver=solver)

    ## Variables
    @variable(m, pflow[1:NLines,TimeChoice:H])
    @variable(m, pgeneration[1:NGenerators,TimeChoice:H])
    @variable(m, storage[1:NNodes,TimeChoice:H] >= 0)
    @variable(m, batterycharge[1:NNodes,TimeChoice:H] >= 0)
    @variable(m, batterydischarge[1:NNodes,TimeChoice:H] >= 0)
    @variable(m, loadshedding[1:NNodes,TimeChoice:H] >= 0)
    @variable(m, productionshedding[1:NNodes,TimeChoice:H] >= 0)
    @variable(m, p_in[TimeChoice:H])
    @variable(m, p_out[TimeChoice:H])

    ## Objective - minimize cost of generation and load shedding
    @objective(m, Min,
        (sum(MargCost[i]*pgeneration[i,u] for i in 1:NGenerators, u = TimeChoice:H)
        + VOLL * sum(loadshedding))
    );

    ## Constraints
    # dynamics
    if TimeChoice == 1
        # if stage 1
        @constraint(m, BatteryDynamics_stage1[n=1:NNodes],
             (storage[n,1] - ini_storage[n]
             - BatteryChargeEfficiency[n] * batterycharge[n,1]
             + batterydischarge[n,1]/BatteryDischargeEfficiency[n]
             == 0)
        )
    else
        # current stage
        @constraint(m, BatteryDynamics_current[n=1:NNodes],
            (storage[n,TimeChoice] - storage_Solution[n,TimeChoice-1,SampleChoice]  # Unique for a sample
             - BatteryChargeEfficiency[n] * batterycharge[n,TimeChoice]
             + batterydischarge[n,TimeChoice]/BatteryDischargeEfficiency[n]
            == 0)
        )
    end
    # future stages
    if TimeChoice < H
        @constraint(m, BatteryDynamics_future[n=1:NNodes, u=TimeChoice+1:H],
            (storage[n,u] - storage[n,u-1]
             - BatteryChargeEfficiency[n] * batterycharge[n,u]
             + batterydischarge[n,u]/BatteryDischargeEfficiency[n]
            == 0)
        )
    end

    # Flow Limits
    @constraint(m, FlowMax[i=1:NLines, u = TimeChoice:H],
        (pflow[i,u] <= SLimit[i])
    )
    @constraint(m, FlowMin[i=1:NLines, u = TimeChoice:H],
        ( - pflow[i,u] <= SLimit[i])
    )

    # Storage Capacity
    @constraint(m, StorageMax[n=1:NNodes, u = TimeChoice:H],
        (storage[n,u] <= BatteryCapacity[n])
    )

    # Charging Capacity
    @constraint(m, BatteryChargeMax[n=1:NNodes, u = TimeChoice:H],
        (batterycharge[n,u] <= BatteryChargeRate[n])
    )

    # Discharging Capacity
    @constraint(m, BatteryDischargeMax[n=1:NNodes, u = TimeChoice:H],
        (batterydischarge[n,u] <= BatteryChargeRate[n])
    )

    # p_in & pflow equality
    @constraint(m, Pin_Flow_equality[u = TimeChoice:H],
        (p_in[u] - pflow[8,u] == 0)
    )
    # p_in & p_out equality
    @constraint(m, Pin_Pout_equality[u = TimeChoice:H],
        (p_in[u] - p_out[u] == 0)
    )

    # Balancing
    # root node
    @constraint(m, Balance1_rootnode[u = TimeChoice:H],
        (sum(pgeneration[g,u] for g = 1:NGenerators)
        + batterydischarge[1,u]+ loadshedding[1,u]
        - productionshedding[1,u] - batterycharge[1,u]
        + sum(pflow[m,u] for m in Children[1])
        == PNetDemand_fix[1,u]
        )
    )
    ## Balancing - usual node
    # @constraint(m, Balance[n=1:NNodes, t in T; n!=0],
    #     (batterydischarge[n,TimeChoice]+ loadshedding[n,TimeChoice]
    #     - productionshedding[n,TimeChoice]- batterycharge[n,TimeChoice]
    #     - pflow[n,TimeChoice]
    #     + sum(pflow[m,TimeChoice] for m in Children[n + 1])
    #     == PNetDemand[n+1,TimeChoice]
    #     )
    # )

    @constraint(m, Balance[n = 1:NNodes, u = TimeChoice:H; n!=0+1 && n!=3+1 && n!=8+1],
        (batterydischarge[n,u]+ loadshedding[n,u]
        - productionshedding[n,u]- batterycharge[n,u]
        - pflow[n-1,u]
        + sum(pflow[m,u] for m in Children[n])
        == PNetDemand_fix[n,u]
        )
    )
    ## Balancing - head node
    @constraint(m, Balance_headnode[n in [8+1], u = TimeChoice:H],
        (batterydischarge[n,u]+ loadshedding[n,u]
        - productionshedding[n,u]- batterycharge[n,u]
        + pflow[n,u]
        + sum(pflow[m,u] for m in Children[n])
        == PNetDemand_fix[n,u]
        )
    )
    ## Balancing - leaf node
    @constraint(m, Balance_leafnode[n in [3+1], u = TimeChoice:H],
        (batterydischarge[n,u]+ loadshedding[n,u]
        - productionshedding[n,u]- batterycharge[n,u]
        - pflow[n,u]
        + sum(pflow[m,u] for m in [4])
        - p_out[u]
        == PNetDemand_fix[n,u]
        )
    )
    ## Generation Limits
    @constraint(m, GenerationMax[i = 1:NGenerators, u = TimeChoice:H],
        (pgeneration[i,u] <= PGenerationMax_fix[i,u])
    )
    @constraint(m, GenerationMin[i = 1:NGenerators, u = TimeChoice:H],
        ( - pgeneration[i,u] <= - PGenerationMin_fix[i,u])
    )
    ## Solve
    # println("Solving problem...")
    status = solve(m)
    # PrintSolution(status)
    # println("Objective value: ", getobjectivevalue(m))

    ## Store Results
    pflow_Solution[:,TimeChoice,SampleChoice] = getvalue(pflow[:,TimeChoice])
    pgeneration_Solution[:,TimeChoice,SampleChoice] = getvalue(pgeneration[:,TimeChoice])
    storage_Solution[:,TimeChoice,SampleChoice] = getvalue(storage[:,TimeChoice])
    batterycharge_Solution[:,TimeChoice,SampleChoice] = getvalue(batterycharge[:,TimeChoice])
    batterydischarge_Solution[:,TimeChoice,SampleChoice] = getvalue(batterydischarge[:,TimeChoice])
    loadshedding_Solution[:,TimeChoice,SampleChoice] = getvalue(loadshedding[:,TimeChoice])
    productionshedding_Solution[:,TimeChoice,SampleChoice] = getvalue(productionshedding[:,TimeChoice])
    p_in_Solution[TimeChoice,SampleChoice] = getvalue(p_in[TimeChoice])
    p_out_Solution[TimeChoice,SampleChoice] = getvalue(p_in[TimeChoice])
    CurrentCost = (sum(MargCost[i]*pgeneration_Solution[i,TimeChoice,SampleChoice] for i in 1:NGenerators)
                    + VOLL * sum(loadshedding_Solution[:,TimeChoice,SampleChoice]))
    return m, CurrentCost
end
## Generate samples scenarios
sample_path = SamplePath(TransProb,NSamples);

## Implementation
SampleCost = Array{Float64}(H,NSamples)
for i = 1:NSamples
    for t = 1:H
        ScenarioChoice = sample_path[:,t,i]
        m, SampleCost[t,i]= CeMPC(t,ScenarioChoice,i);
        println("   cost of stage ",t,", sample ",i , SampleCost[t,i])
    end
    println("Total cost of sample ", i, " : ", sum(SampleCost[:,i]))
    # print(sum(CurrentCost[:,i]))
end
