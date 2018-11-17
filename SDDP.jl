using JuMP, Gurobi

# Here we load the SDDP data
include("LoadDataSDDP.jl")

# The struct to Store the solutions

struct Solutions
    " Struct that stores the decision variables"
    # Saved as Arrays
    pflow::Dict{Int64, Float64}
    storage::Dict{Int64, Float64}
    batterycharge::Dict{Int64, Float64}
    batterydischarge::Dict{Int64, Float64}
    loadshedding::Dict{Int64, Float64}
    productionshedding::Dict{Int64, Float64}
    pgeneration

    " Struct that stores the dual multipliers"
    # saved as Dictionaries
    BatteryDynamics::Dict{Int64, Float64}
    Balance::Dict{Int64, Float64}
    Pin_fix::Dict{Int64, Float64}
    Pout_fix::Dict{Array{Int64,1}, Float64}
    Pin_Flow_equality::Dict{Int64, Float64}
    FlowMax::Dict{Int64, Float64}
    FlowMin::Dict{Int64, Float64}
    StorageMax::Dict{Int64, Float64}
    BatteryChargeMax::Dict{Int64, Float64}
    BatteryDischargeMax::Dict{Int64, Float64}
    Cuts
    GenerationMax
    GenerationMin

    "Optimal Value"
    OptimalValue::Float64
end

# Function that translate JuMP.Dict type of
# the dual multipliers to Dict type Dict{Int64,Float64}
function CreateDictionaryV(ArrayChoice)
    Num = length(ArrayChoice);
    Dictionary = Dict{Int64,Float64}();
    for i = 1:Num
        for n in collect(keys(ArrayChoice[i]))
            Dictionary[n[1]] = ArrayChoice[i][n[1]];
        end
    end
    return Dictionary
end

# Function that translate JuMP.Dict type of
# the dual multipliers to Dict type Dict{Array{Int64,1},Float64}

function CreateDictionaryA(ArrayChoice)
    Num = length(ArrayChoice);
    Dictionary = Dict{Array{Int64,1},Float64}();
    for i = 1:Num
        for n in collect(keys(ArrayChoice[i]))
            Dictionary[[n[1],n[2]]] = ArrayChoice[i][n[1],n[2]];
        end
    end
    return Dictionary
end

# Function that determines wether the obtained cut is repeated
# 0 is not repeated, 1 yes
function UsefulCut(LayerChoice, TimeChoice, etemp, Etemp)
    i = 1;
    boo = 0;
    while boo == 0 && i <= length(eCut[LayerChoice, TimeChoice])
        Difference = 0.0;
        for n in LayerNodes[LayerChoice]
            Difference = Difference + abs(ECut[LayerChoice, TimeChoice, n][i] - Etemp[n]);
        end
        if abs(eCut[LayerChoice, TimeChoice][i] - etemp) < 10e-6 &&  Difference < 10e-6
            boo = 1;
        else
            i = i+1;
        end
    end
    return boo
end

# The NLDS algorithm
function NLDS(LayerChoice, SampleChoice, TimeChoice, OutcomeChoice, iter)
    " NLDS(t,k) "
    " LayerChoice -- Refers to the current Layer"
    " SampleChoice -- Refers to the index of the sample of the Monte Carlo"
    " TimeChoice -- Refers to the time stage t"
    " OutcomeChoice -- Refers to the outcome k of time stage t"

    # Definition of the model
    model = Model(solver=GurobiSolver())
    # Variables
    @variable(model, pflow[LayerLines[LayerChoice]]  );
    @variable(model, storage[LayerNodes[LayerChoice]] >= 0);
    @variable(model, batterycharge[LayerNodes[LayerChoice]] >= 0);
    @variable(model, batterydischarge[LayerNodes[LayerChoice]] >= 0);
    @variable(model, loadshedding[LayerNodes[LayerChoice]] >= 0 );
    @variable(model, productionshedding[LayerNodes[LayerChoice]] >= 0);
    @variable(model, p_in[HeadNodes[LayerChoice]]);
    @variable(model, p_out[n in LeafNodes[LayerChoice], LeafChildren[LayerChoice,n]] );
    @variable(model, theta >= 0 );

    # Definition of the objective function and additional variables, constraints
    # in case we are in the root Node

    if LayerChoice == 1
        @variable(model, pgeneration[1:NGenerators] );
        @objective(model, Min,  sum(MargCost[g]*pgeneration[g] for g=1:NGenerators) + VOLL*sum(loadshedding[n] for n in LayerNodes[LayerChoice]) + theta)

        # Balancing - root node
        @constraint(model, Balance_rootnode, sum(pgeneration[g] for g = 1:NGenerators) + batterydischarge[1] + loadshedding[1] - productionshedding[1] - batterycharge[1] + sum(pflow[m] for m in Children[1]) ==  PNetDemand[1,TimeChoice][OutcomeChoice] );

        # Generation Limits
        @constraint(model, GenerationMax[g = 1:NGenerators], pgeneration[g] <= PGenerationMax[g,TimeChoice][OutcomeChoice] );
        @constraint(model, GenerationMin[g = 1:NGenerators],  -pgeneration[g] <= -PGenerationMin[g,TimeChoice][OutcomeChoice]);

    else
        @objective(model, Min,  VOLL*sum(loadshedding[n] for n in LayerNodes[LayerChoice]) + theta)
    end

    ##### Constraints

    # Battery Constraints
    if TimeChoice == 1
        @constraint(model,  BatteryDynamics[n in LayerNodes[LayerChoice]], storage[n] - BatteryChargeEfficiency[n] * batterycharge[n] + batterydischarge[n]/BatteryDischargeEfficiency[n] - ini_storage[n] == 0 )
    else
        @constraint(model,  BatteryDynamics[n in LayerNodes[LayerChoice]], storage[n] - BatteryChargeEfficiency[n] * batterycharge[n] + batterydischarge[n]/BatteryDischargeEfficiency[n] - storageTrials[LayerChoice, n, TimeChoice-1, SampleChoice, iter] == 0 )
    end

    # Balance Constraints
    # Balancing - usual nodes
    @constraint(model, Balance[n in setdiff(LayerNodes[LayerChoice], union( [1], HeadNodes[LayerChoice], LeafNodes[LayerChoice] )) ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m]  for m in Children[n]) - pflow[n-1]  ==  PNetDemand[n,TimeChoice][OutcomeChoice] )  );

    # Balancing - head node
    @constraint(model, Balance_headnode[n in setdiff(HeadNodes[LayerChoice], LeafNodes[LayerChoice]) ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m] for m in Children[n]) + pflow[n-1] ==  PNetDemand[n,TimeChoice][OutcomeChoice] ) );

    # Balancing - leaf node
    @constraint(model, Balance_leafnode[n in setdiff(LeafNodes[LayerChoice], HeadNodes[LayerChoice]) ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m] for m in Children[n] if (issubset([m+1], LeafChildren[LayerChoice,n]) == false) )
    - pflow[n-1] - sum(p_out[n,j] for j in LeafChildren[LayerChoice,n]) ==  PNetDemand[n,TimeChoice][OutcomeChoice] )  );

    # Balancing - head-leaf node
    @constraint(model, Balance_headleafnode[n in intersect(LeafNodes[LayerChoice], HeadNodes[LayerChoice]) ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m] for m in Children[n] if issubset([m+1], LeafChildren[LayerChoice,n]) == false )
    + pflow[n-1] - sum(p_out[n,j] for j in LeafChildren[LayerChoice,n]) ==  PNetDemand[n,TimeChoice][OutcomeChoice] ) );

    # p_in constraint
    @constraint(model, Pin_fix[n in HeadNodes[LayerChoice]], p_in[n] == p_in_data[TimeChoice,OutcomeChoice] );

    # p_out constraint
    @constraint(model, Pout_fix[n in LeafNodes[LayerChoice], m in LeafChildren[LayerChoice,n]], p_out[n,m] == p_out_data[TimeChoice,OutcomeChoice] );

    # p_in & pflow equality constraint
    @constraint(model, Pin_Flow_equality[n in HeadNodes[LayerChoice]], p_in[n] - pflow[n-1] == 0);

    # Max Flow Limit contraint
    @constraint(model, FlowMax[n in LayerLines[LayerChoice]], pflow[n] <= SLimit[n]  );

    # Min Flow Limit contraint
    @constraint(model, FlowMin[n in LayerLines[LayerChoice]], -pflow[n] <= SLimit[n]  );

    # Storage Capacity constraint
    @constraint(model, StorageMax[n in LayerNodes[LayerChoice] ], storage[n] <= BatteryCapacity[n] );

    # Charging Capacity constraint
    @constraint(model, BatteryChargeMax[n in LayerNodes[LayerChoice] ], batterycharge[n] <= BatteryChargeRate[n]);

    # Discharging Capacity constraint
    @constraint(model, BatteryDischargeMax[n in LayerNodes[LayerChoice]], batterydischarge[n] <= BatteryChargeRate[n]);

    # Optimality cuts and theta = zero in case TimeChoice is H
    if TimeChoice == H
        @constraint(model, theta == 0 )
    else
        @constraint(model, Cuts[i = 1:length(eCut[LayerChoice, TimeChoice])],
        theta  >= eCut[LayerChoice, TimeChoice][i] - sum(ECut[LayerChoice, TimeChoice, n][i]*storage[n] for n in LayerNodes[LayerChoice]) );
    end

    ##### Solve the Model
    solve(model)

    ##### Here we return the results
    # Here we translate JuMP.Dict type of dual Multipliers & optimal decisions to Dict Type
    BatteryDynamicsMultipliers =  CreateDictionaryV( push!([], getdual(BatteryDynamics)) );
    BalanceMultipliers = CreateDictionaryV( push!([], getdual(Balance), getdual(Balance_headnode), getdual(Balance_leafnode), getdual(Balance_headleafnode) ));
    Pin_fixMultipliers = CreateDictionaryV( push!([], getdual(Pin_fix)) );
    Pout_fixMultipliers = CreateDictionaryA( push!([], getdual(Pout_fix)) );
    FlowMaxMultipliers = CreateDictionaryV( push!([], getdual(FlowMax)) );
    FlowMinMultipliers = CreateDictionaryV( push!([], getdual(FlowMin)) );
    Pin_Flow_equalityMultipliers = CreateDictionaryV( push!([], getdual(Pin_Flow_equality)) );
    StorageMaxMultipliers = CreateDictionaryV( push!([], getdual(StorageMax)) );
    BatteryChargeMaxMultipliers = CreateDictionaryV( push!([], getdual(BatteryChargeMax)) );
    BatteryDischargeMaxMultipliers = CreateDictionaryV( push!([], getdual(BatteryDischargeMax)) );

    # Return certain data depending on the LayerChoice and TimeChoice
    if LayerChoice == 1
        BalanceMultipliers[1] = getdual(Balance_rootnode);
        GenerationMaxMultipliers = CreateDictionaryV( push!([], getdual(GenerationMax)) );
        GenerationMinMultipliers = CreateDictionaryV( push!([], getdual(GenerationMin)) );
        if TimeChoice == H
            return Solutions(  CreateDictionaryV(push!([], getvalue(pflow))), CreateDictionaryV(push!([],getvalue(storage))),
            CreateDictionaryV(push!([],getvalue(batterycharge))), CreateDictionaryV(push!([],getvalue(batterydischarge))),
            CreateDictionaryV(push!([],getvalue(loadshedding))), CreateDictionaryV( push!([],getvalue(productionshedding))),
            getvalue(pgeneration)[:], BatteryDynamicsMultipliers,
            BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers, Pin_Flow_equalityMultipliers, FlowMaxMultipliers,
            FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers, BatteryDischargeMaxMultipliers, "NoCuts", GenerationMaxMultipliers,
            GenerationMinMultipliers, getobjectivevalue(model) )
        else
            return Solutions(CreateDictionaryV(push!([], getvalue(pflow))), CreateDictionaryV(push!([],getvalue(storage))),
            CreateDictionaryV(push!([],getvalue(batterycharge))), CreateDictionaryV(push!([],getvalue(batterydischarge))),
            CreateDictionaryV(push!([],getvalue(loadshedding))), CreateDictionaryV( push!([],getvalue(productionshedding))),
            getvalue(pgeneration)[:], BatteryDynamicsMultipliers,
            BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers, Pin_Flow_equalityMultipliers, FlowMaxMultipliers,
            FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers, BatteryDischargeMaxMultipliers, getdual(Cuts), GenerationMaxMultipliers,
            GenerationMinMultipliers, getobjectivevalue(model) )
        end
    else
        if TimeChoice == H
            return Solutions(CreateDictionaryV(push!([], getvalue(pflow))), CreateDictionaryV(push!([],getvalue(storage))),
            CreateDictionaryV(push!([],getvalue(batterycharge))), CreateDictionaryV(push!([],getvalue(batterydischarge))),
            CreateDictionaryV(push!([],getvalue(loadshedding))), CreateDictionaryV( push!([],getvalue(productionshedding))),
            "No pgeneration", BatteryDynamicsMultipliers, BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers,
            Pin_Flow_equalityMultipliers, FlowMaxMultipliers, FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers,
            BatteryDischargeMaxMultipliers, "NoCuts", "No pgeneration", "No pgeneration", getobjectivevalue(model) )
        else
            return Solutions(CreateDictionaryV(push!([], getvalue(pflow))), CreateDictionaryV(push!([],getvalue(storage))),
            CreateDictionaryV(push!([],getvalue(batterycharge))), CreateDictionaryV(push!([],getvalue(batterydischarge))),
            CreateDictionaryV(push!([],getvalue(loadshedding))), CreateDictionaryV( push!([],getvalue(productionshedding))),
            "No pgeneration", BatteryDynamicsMultipliers, BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers,
            Pin_Flow_equalityMultipliers, FlowMaxMultipliers, FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers,
            BatteryDischargeMaxMultipliers, getdual(Cuts), "No pgeneration", "No pgeneration", getobjectivevalue(model) )
        end
    end
end


# Function for the forward pass
function ForwardPass(K, LayerChoice, iter)
    " K -- Refers to the number of Monte Carlo Samples"
    " LayerChoice -- Refers to the current Layer"
    " iter -- Refers to the current iteration"

    # Solve the First Stage
    FirstNLDS = NLDS(LayerChoice, 1, 1, 1, iter);

    # Store The Lower Bound
    LowerBound[LayerChoice, iter] = FirstNLDS.OptimalValue
    # Store the results
    for n in LayerLines[LayerChoice]
        pflowTrials[LayerChoice, n, 1, :, iter] = FirstNLDS.pflow[n]*ones(K);
    end
    for n in LayerNodes[LayerChoice]
        storageTrials[LayerChoice, n, 1, :, iter] = FirstNLDS.storage[n]*ones(K);
        batterychargeTrials[LayerChoice, n, 1, :, iter] = FirstNLDS.batterycharge[n]*ones(K);
        batterydischargeTrials[LayerChoice, n, 1, :, iter] = FirstNLDS.batterydischarge[n]*ones(K);
        loadsheddingTrials[LayerChoice, n, 1, : , iter] = FirstNLDS.loadshedding[n]*ones(K);
        productionsheddingTrials[LayerChoice, n, 1, :, iter] = FirstNLDS.productionshedding[n]*ones(K);
    end
    if LayerChoice == 1
        pgenerationTrials[1:NGenerators, 1, :, iter] = repmat(FirstNLDS.pgeneration,1,K)
    end

    # Solve the remaining stages
    for SampleChoice = 1:K
        for TimeChoice = 2:H
            # Solve the Stage
            path[TimeChoice,SampleChoice,iter] = rand(1:NLattice[TimeChoice]);
            NestedLDS = NLDS(LayerChoice, SampleChoice, TimeChoice, path[TimeChoice, SampleChoice, iter], iter);

            # Store the results
            for n in LayerLines[LayerChoice]
                pflowTrials[LayerChoice, n, TimeChoice, SampleChoice, iter] = NestedLDS.pflow[n];
            end
            for n in LayerNodes[LayerChoice]
                batterychargeTrials[LayerChoice, n, TimeChoice, SampleChoice, iter] = NestedLDS.batterycharge[n];
                storageTrials[LayerChoice, n, TimeChoice, SampleChoice, iter] = NestedLDS.storage[n] ;
                batterydischargeTrials[LayerChoice, n, TimeChoice, SampleChoice, iter] = NestedLDS.batterydischarge[n];
                loadsheddingTrials[LayerChoice, n, TimeChoice, SampleChoice, iter] = NestedLDS.loadshedding[n] ;
                productionsheddingTrials[LayerChoice, n, TimeChoice, SampleChoice, iter] = NestedLDS.productionshedding[n];
            end
            if LayerChoice == 1
                pgenerationTrials[1:NGenerators, TimeChoice, SampleChoice, iter] = NestedLDS.pgeneration ;
            end
        end
    end
end

# Function for the BackwardPass
function BackwardPass(K, LayerChoice, iter)
    " K -- Refers to the number of Monte Carlo Samples"
    " LayerChoice -- Refers to the current Layer"
    " iter -- Refers to the current iteration"
    for t = 1:H-1
        for SampleChoice = 1:K
            etemp = 0.0;
            Etemp = Dict{Int64,Float64}(n => 0.0 for n in LayerNodes[LayerChoice]);
            for OutcomeChoice = 1:NLattice[H+1-t]
                # The NLDS is solved
                NesLDS = NLDS(LayerChoice, SampleChoice, H-t+1, OutcomeChoice, iter);
                # The E coefficients are calculated
                # The e coeficients are calculated
                for n in HeadNodes[LayerChoice]
                    etemp = etemp + NesLDS.Pin_fix[n]*p_in_data[H-t+1,OutcomeChoice] ;
                end
                for n in LayerNodes[LayerChoice]
                    etemp = etemp + NesLDS.Balance[n]*PNetDemand[n,H-t+1][OutcomeChoice] + NesLDS.StorageMax[n]*BatteryCapacity[n] + NesLDS.BatteryChargeMax[n]*BatteryChargeRate[n] + NesLDS.BatteryDischargeMax[n]*BatteryChargeRate[n];
                    Etemp[n] = Etemp[n] - TransProb[LayerChoice,H-t][OutcomeChoice]*NesLDS.BatteryDynamics[n];
                end
                for n in LayerLines[LayerChoice]
                    etemp = etemp + NesLDS.FlowMax[n]*SLimit[n] + 2*NesLDS.FlowMin[n]*SLimit[n] ;
                end
                for n in LeafNodes[LayerChoice], m in LeafChildren[LayerChoice,n]
                    etemp = etemp + NesLDS.Pout_fix[[n, m]]*p_out_data[H-t+1,OutcomeChoice];
                end
                if LayerChoice == 1
                    etemp = etemp + sum( NesLDS.GenerationMax[g]*PGenerationMax[g, H-t+1][OutcomeChoice] for g = 1:NGenerators ) - sum( NesLDS.GenerationMin[g]*PGenerationMin[g,H-t+1][OutcomeChoice] for g = 1:NGenerators )
                end
                if t > 1
                    etemp = etemp + sum( NesLDS.Cuts[:].*eCut[LayerChoice, H+1-t][:]) ;
                end
                etemp = TransProb[LayerChoice,H-t][OutcomeChoice]*etemp;
            end
            # Here we check if the cut is useful
            if UsefulCut(LayerChoice, H-t, etemp, Etemp) == 0
                for n in LayerNodes[LayerChoice]
                    push!(ECut[LayerChoice, H-t, n], Etemp[n]);
                end
                push!(eCut[LayerChoice, H-t], etemp);
            end
        end
    end
end
