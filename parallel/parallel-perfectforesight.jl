"""
Define the function of perfect foresight with JuMP for parallel computing
"""

## define solver
solver = GurobiSolver(LogToConsole=0#=, LogFile="log/parallel-perfectforesight.log"=#)

function PerfectForesight(RealPath, solutions, idx=NaN)
    "
    implement perfect foresight policy
    Args:
        - RealPath: Real path of the 'current' sample
        - solutions: an instance of the struct Solutions
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
    ## quadratic optimization
    # @objective(m, Min,
    #     (sum(MargCost[i]*(pgeneration[i,u])^2 for i = 1:NGenerators, u = 1:H)
    #     + VOLL * sum(loadshedding[n,u]^2 for n=1:NNodes,u=1:H))
    # );

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
    # @time status = solve(m);
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
    if isnan(idx)
        println("Total cost is $(sum(solutions.StageCost[:]))")
    else
        println("Sample $idx, Total cost is $(sum(solutions.StageCost[:]))")
    end
    return #solutions.StageCost[:]
end

function PerfectForesight_all(idx)
    "
    Implement PerfectForesight over a sample
    "
    solutions = Solutions()  # define solution struct
    RealPath = SamplePath(TransProb);  # generate a sample
    tic();
    PerfectForesight(RealPath, solutions, idx)
    solutions.IterationTime[1] = toc();  # store the time in the first element
    return solutions
end
