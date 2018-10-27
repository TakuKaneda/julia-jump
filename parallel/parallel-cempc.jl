"""
Define the function of ce-mpc with JuMP for parallel computing
"""

## define solver
solver = GurobiSolver(LogToConsole=0, LogFile="log/parallel-cempc.log")

function ComputeExpectedParameters(TimeChoice,ScenarioChoice)
    "
    Compute expected value of stochastic parameters of the problem, i.e.
    * PNetDemand_fix
    * PGenerationMax_fix
    * PGenerationMin_fix
    Args:
        - TimeChoice: 'current' stage
        - ScenarioChoice[l]: vector of outcomes for layer l
    **NOTE: return the current stage cost **
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
function CeMPC(TimeChoice, RealPath, solutions)

    "
    implement ce-mpc (certainty-equivalent MPC)
    Args:
        - TimeChoice: 'current' stage
        - RealPath: Real path of the 'current' sample
        - solutions: an instance of the struct Solutions
    note: results will be stored in `solutions` thus return is Void
    "
    ## Compute expected value of stochastic parameters
    ScenarioChoice = RealPath[:,TimeChoice]
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
    # @variable(m, p_in[TimeChoice:H])
    # @variable(m, p_out[TimeChoice:H])

    ## Objective - minimize cost of generation and load shedding
    @objective(m, Min,
        (sum(MargCost[i]*pgeneration[i,u] for i = 1:NGenerators, u = TimeChoice:H)
        + VOLL * sum(loadshedding))
    );
    ## quadratic optimization
    # @objective(m, Min,
    #     (sum(MargCost[i]*(pgeneration[i,u])^2 for i = 1:NGenerators, u = TimeChoice:H)
    #     + VOLL * sum(loadshedding[n,u]^2 for n=1:NNodes,u=TimeChoice:H))
    # );

    ## Constraints
    # dynamics
    if TimeChoice == 1
        # if stage 1
        @constraint(m, BatteryDynamics_stage1[n=1:NNodes],
             (storage[n,1] - ini_storage[n]
             - BatteryChargeEfficiency[n] * batterycharge[n,1]
             + batterydischarge[n,1]/BatteryDischargeEfficiency[n]
             == 0)
        );
    else
        # current stage
        @constraint(m, BatteryDynamics_current[n=1:NNodes],
            (storage[n,TimeChoice] - solutions.storage[n,TimeChoice-1]  # Unique for a sample
             - BatteryChargeEfficiency[n] * batterycharge[n,TimeChoice]
             + batterydischarge[n,TimeChoice]/BatteryDischargeEfficiency[n]
            == 0)
        );
    end
    # future stages
    if TimeChoice < H
        @constraint(m, BatteryDynamics_future[n=1:NNodes, u=TimeChoice+1:H],
            (storage[n,u] - storage[n,u-1]
             - BatteryChargeEfficiency[n] * batterycharge[n,u]
             + batterydischarge[n,u]/BatteryDischargeEfficiency[n]
            == 0)
        );
    end

    # Flow Limits
    @constraint(m, FlowMax[i=1:NLines, u = TimeChoice:H],
        (pflow[i,u] <= SLimit[i])
    );
    @constraint(m, FlowMin[i=1:NLines, u = TimeChoice:H],
        ( - pflow[i,u] <= SLimit[i])
    );

    # Storage Capacity
    @constraint(m, StorageMax[n=1:NNodes, u = TimeChoice:H],
        (storage[n,u] <= BatteryCapacity[n])
    );

    # Charging Capacity
    @constraint(m, BatteryChargeMax[n=1:NNodes, u = TimeChoice:H],
        (batterycharge[n,u] <= BatteryChargeRate[n])
    );

    # Discharging Capacity
    @constraint(m, BatteryDischargeMax[n=1:NNodes, u = TimeChoice:H],
        (batterydischarge[n,u] <= BatteryChargeRate[n])
    );
    #
    # # p_in & pflow equality
    # @constraint(m, Pin_Flow_equality[u = TimeChoice:H],
    #     (p_in[u] - pflow[8,u] == 0)
    # );
    # # p_in & p_out equality
    # @constraint(m, Pin_Pout_equality[u = TimeChoice:H],
    #     (p_in[u] - p_out[u] == 0)
    # );

    # Balancing
    # root node
    @constraint(m, Balance1_rootnode[u = TimeChoice:H],
        (sum(pgeneration[g,u] for g = 1:NGenerators)
        + batterydischarge[1,u]+ loadshedding[1,u]
        - productionshedding[1,u] - batterycharge[1,u]
        + sum(pflow[m,u] for m in Children[1])
        == PNetDemand_fix[1,u]
        )
    );
    # Balancing - usual nodes
    @constraint(m, Balance[n = 1:NNodes, u = TimeChoice:H; n!=0+1 #=&& n!=3+1 && n!=8+1=#],
        (batterydischarge[n,u]+ loadshedding[n,u]
        - productionshedding[n,u]- batterycharge[n,u]
        - pflow[n-1,u]
        + sum(pflow[m,u] for m in Children[n])
        == PNetDemand_fix[n,u]
        )
    );
    # Balancing - head node
    # @constraint(m, Balance_headnode[n in [8+1], u = TimeChoice:H],
    #     (batterydischarge[n,u]+ loadshedding[n,u]
    #     - productionshedding[n,u]- batterycharge[n,u]
    #     + pflow[n-1,u]
    #     + sum(pflow[m,u] for m in Children[n])
    #     == PNetDemand_fix[n,u]
    #     )
    # );
    # # Balancing - leaf node
    # @constraint(m, Balance_leafnode[n in [3+1], u = TimeChoice:H],
    #     (batterydischarge[n,u]+ loadshedding[n,u]
    #     - productionshedding[n,u]- batterycharge[n,u]
    #     - pflow[n-1,u]
    #     + sum(pflow[m,u] for m in [4])
    #     - p_out[u]
    #     == PNetDemand_fix[n,u]
    #     )
    # );
    # Generation Limits
    @constraint(m, GenerationMax[i = 1:NGenerators, u = TimeChoice:H],
        (pgeneration[i,u] <= PGenerationMax_fix[i,u])
    );
    @constraint(m, GenerationMin[i = 1:NGenerators, u = TimeChoice:H],
        ( - pgeneration[i,u] <= - PGenerationMin_fix[i,u])
    );

    ## Solve
    @time status = solve(m);

    ## Store Results
    solutions.pflow[:,TimeChoice] = getvalue(pflow[:,TimeChoice])
    solutions.pgeneration[:,TimeChoice] = getvalue(pgeneration[:,TimeChoice])
    solutions.storage[:,TimeChoice] = getvalue(storage[:,TimeChoice])
    solutions.batterycharge[:,TimeChoice] = getvalue(batterycharge[:,TimeChoice])
    solutions.batterydischarge[:,TimeChoice] = getvalue(batterydischarge[:,TimeChoice])
    solutions.loadshedding[:,TimeChoice] = getvalue(loadshedding[:,TimeChoice])
    solutions.productionshedding[:,TimeChoice] = getvalue(productionshedding[:,TimeChoice])
    # solutions.p_in[TimeChoice] = getvalue(p_in[TimeChoice])
    # solutions.p_out[TimeChoice] = getvalue(p_out[TimeChoice])
    solutions.StageCost[TimeChoice] = (sum(MargCost[i]*solutions.pgeneration[i,TimeChoice] for i in 1:NGenerators)
                    + VOLL * sum(solutions.loadshedding[:,TimeChoice]))
    return solutions.StageCost[TimeChoice]  # return the current stage cost
end
