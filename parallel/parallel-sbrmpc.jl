"""
Define the function of sbr-mpc with JuMP for parallel computing
"""

## define solver
solver = GurobiSolver(LogToConsole=0, LogFile="log/parallel-sbrmpc.log")

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


function SbrMPC(TimeChoice, RealPath, solutions, idx=NaN)
    "
    implement sbr-mpc (scenario-based robust MPC)
    Args:
        - TimeChoice: 'current' stage
        - RealPath: Real path of the 'current' sample
        - solutions: an instance of the struct Solutions
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
    # @time status = solve(m);
    status = solve(m);

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
    if isnan(idx)
        println("Stage $TimeChoice, StageCost is $(solutions.StageCost[TimeChoice])")
    else
        println("Sample $idx, Stage $TimeChoice, StageCost is $(solutions.StageCost[TimeChoice])")
    end
    return #solutions.StageCost[TimeChoice]  # return the current stage cost
end

function SbrMPC_all(idx)
    "
    Implement sbr MPC over a sample
    "
    solutions = Solutions()  # define solution struct
    RealPath = SamplePath(TransProb);  # generate a sample
    println("Sample $idx Start")
    for  t = 1:H
        tic();
        SbrMPC(t, RealPath, solutions, idx)
        solutions.IterationTime[t] = toc();
    end
    println("Sample $idx Finish")
    return solutions
end
