### Implement ce-mpc with JuMP

using DataFrames, JuMP, Gurobi, CSV, JSON
solver = GurobiSolver()

## Read CSV
lines_df = CSV.read("data/twolayer-lines.csv")
nodes_df = CSV.read("data/twolayer-nodes.csv")
generators_df = CSV.read("data/twolayer-generators.csv")
netdemand_df = CSV.read("data/sample-netdemand.csv")

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
# print(netdemand_df)

## Problem Parameters
# generators
Generators = generators_df[:GeneratorID]
MargCost = generators_df[:MargCost]
PMax = generators_df[:PMax]
PMin = generators_df[:PMin]
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
# net demand
NetDemand = Array(netdemand_df)
# other parameters
NLines = nrow(lines_df)
NNodes = nrow(nodes_df)
NGenerators = nrow(generators_df)
VOLL = 5000
H = 24
T = 1:H

## implement certainty-equivalent Model Predictive Control
function CeMPC(t)
    ## Build model
    m = Model(solver=solver)

    ## Variables
    @variable(m, pflow[Lines,t:H])
    @variable(m, pgeneration[Generators,t:H])
    @variable(m, storage[Nodes,t:H] >= 0)
    @variable(m, batterycharge[Nodes,t:H] >= 0)
    @variable(m, batterydischarge[Nodes,t:H] >= 0)
    @variable(m, loadshedding[Nodes,t:H] >= 0)
    @variable(m, productionshedding[Nodes,t:H] >= 0)
    @variable(m, p_in[t:H])
    @variable(m, p_out[t:H])

    ## Objective - minimize cost of generation and load shedding
    @objective(m, Min,
        (sum(MargCost[i]*pgeneration[Generators[i],t] for i in 1:NGenerators, t in T)
        + VOLL * sum(loadshedding))
    )

    ## Constraints
    # dynamics
    @constraint(m, BatteryDynamics_t1[n in Nodes],
         (storage[n,1] - ini_storage[n+1]
         - BatteryChargeEfficiency[n+1] * batterycharge[n,1]
         + batterydischarge[n,1]/BatteryDischargeEfficiency[n+1]
         == 0)
    )

    @constraint(m, BatteryDynamics[n in Nodes, t=2:H],
        (storage[n,t] - storage[n,t-1]
         - BatteryChargeEfficiency[n+1] * batterycharge[n,t]
         + batterydischarge[n,t]/BatteryDischargeEfficiency[n+1]
        == 0)
    )
    ## Flow Limits
    @constraint(m, FlowMax[i in Lines, t in T],
        (pflow[i,t] <= SLimit[i])
    )
    @constraint(m, FlowMin[i in Lines, t in T],
        ( - pflow[i,t] <= SLimit[i])
    )

    ## Generation Limits
    @constraint(m, GenerationMax[i in 1:NGenerators, t in T],
        (pgeneration[Generators[i],t] <= PMax[i])
    )
    @constraint(m, GenerationMin[i in 1:NGenerators, t in T],
        ( - pgeneration[Generators[i],t] <= - PMin[i])
    )

    ## Storage Capacity
    @constraint(m, StorageMax[n in Nodes, t in T],
        (storage[n,t] <= BatteryCapacity[n+1])
    )
    ## Charging Capacity
    @constraint(m, BatteryChargeMax[n in Nodes, t in T],
        (batterycharge[n,t] <= BatteryChargeRate[n+1])
    )
    ## Discharging Capacity
    @constraint(m, BatteryDischargeMax[n in Nodes, t in T],
        (batterydischarge[n,t] <= BatteryChargeRate[n+1])
    )
    # p_in & pflow equality
    @constraint(m, Pin_Flow_equality[t in T],
        (p_in[t] - pflow[8,t] == 0)
    )
    # p_in & p_out equality
    @constraint(m, Pin_Pout_equality[t in T],
        (p_in[t] - p_out[t] == 0)
    )

    ## Balancing - root node
    @constraint(m, Balance1_node0[t in T],
        (sum(pgeneration[g,t] for g in Generators)
        + batterydischarge[0,t]+ loadshedding[0,t]
        - productionshedding[0,t] - batterycharge[0,t]
        + sum(pflow[m,t] for m in Children[0+1])
        == NetDemand[0+1,t]
        )
    )
    ## Balancing - usual node
    # @constraint(m, Balance[n in Nodes, t in T; n!=0],
    #     (batterydischarge[n,t]+ loadshedding[n,t]
    #     - productionshedding[n,t]- batterycharge[n,t]
    #     - pflow[n,t]
    #     + sum(pflow[m,t] for m in Children[n + 1])
    #     == NetDemand[n+1,t]
    #     )
    # )

    @constraint(m, Balance[n in Nodes, t in T; n!=0 && n!=3 && n!=8],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        - pflow[n,t]
        + sum(pflow[m,t] for m in Children[n + 1])
        == NetDemand[n+1,t]
        )
    )
    ## Balancing - head node
    @constraint(m, Balance_headnode[n in [8], t in T],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        + pflow[n,t]
        + sum(pflow[m,t] for m in Children[n + 1])
        == NetDemand[n+1,t]
        )
    )
    ## Balancing - leaf node
    @constraint(m, Balance_leafnode[n in [3], t in T],
        (batterydischarge[n,t]+ loadshedding[n,t]
        - productionshedding[n,t]- batterycharge[n,t]
        - pflow[n,t]
        + sum(pflow[m,t] for m in [4])
        - p_out[t]
        == NetDemand[n+1,t]
        )
    )

    ## Solve
    println("Solving problem...")
    status = solve(m)
    # PrintSolution(status)
    println("Objective value: ", getobjectivevalue(m))
    return m
end

m = CeMPC()
