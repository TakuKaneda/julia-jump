using JuMP, Gurobi

# Here we load the SDDP data
include("LoadDataSDDP.jl")

### TAKU add
# define the cuts
# eCut[l, t][k] has a vector of cut-RHS (`e^l_{t,k}` of NLDS of layer l at stage t, outcome k)
# ECut[n, t][k] has a vector of cut-coefficients of node n (`E^l_{n,t,k}` of NLDS of layer l at stage t, outcome k)
#   note that there are no cuts at stage H
eCut = Array{Array}(NLayers,H-1)
ECut = Array{Array}(NNodes,H-1)
for t = 1:H-1
    for l = 1:NLayers
        eCut[l,t] = Array{Any}(NLattice[t]) # number of lattice
    end
    for n = 1:NNodes
        ECut[n,t] = Array{Any}(NLattice[t]) # number of lattice
    end
end
### end TAKU add

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
# function UsefulCut(LayerChoice, TimeChoice, etemp, Etemp)
#     i = 1;
#     boo = 0;
#     while boo == 0 && i <= length(eCut[LayerChoice, TimeChoice])
#         Difference = 0.0;
#         for n in LayerNodes[LayerChoice]
#             Difference = Difference + abs(ECut[LayerChoice, TimeChoice, n][i] - Etemp[n]);
#         end
#         if abs(eCut[LayerChoice, TimeChoice][i] - etemp) < 10e-6 &&  Difference < 10e-6
#             boo = 1;
#         else
#             i = i+1;
#         end
#     end
#     return boo
# end

### TAKU add
# Function that determines wether the obtained cut is repeated
function UsefulCut(LayerChoice, TimeChoice, OutcomeChoice, etemp, Etemp)
    "Function that determines wether the obtained cut is repeated
    return  `true` : the cut is useful
            `false`: otherwise
    "
    i = 1;
    boo = true;
    while boo == true && i <= length(eCut[LayerChoice, TimeChoice][OutcomeChoice])
        Difference = sum(abs(ECut[n, TimeChoice][OutcomeChoice][i] - Etemp[n]) for n in LayerNodes[LayerChoice]) # compare E
        Difference += abs(eCut[LayerChoice, TimeChoice][OutcomeChoice][i] - etemp) # compare e
        if Difference < 10e-6
            boo = false; # the cut is repeated i.e. unuseful
        else
            i = i+1; # go to the next cut that will be compared
        end
    end
    return boo
end
### end TAKU add

# The NLDS algorithm
function NLDS(LayerChoice, SampleChoice, TimeChoice, OutcomeChoice, iter)
    " NLDS(t,k) "
    " LayerChoice -- Refers to the current Layer"
    " SampleChoice -- Refers to the index of the sample of the Monte Carlo"
    " TimeChoice -- Refers to the time stage t"
    " OutcomeChoice -- Refers to the outcome k of time stage t"

    # Definition of the model
    model = Model(solver=GurobiSolver(LogToConsole=0))
    # model = Model(solver=GurobiSolver())
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
        @constraint(model,  BatteryDynamics[n in LayerNodes[LayerChoice]], storage[n] - BatteryChargeEfficiency[n] * batterycharge[n] + batterydischarge[n]/BatteryDischargeEfficiency[n] - ini_storage[n] == 0.0 )
    else
        @constraint(model,  BatteryDynamics[n in LayerNodes[LayerChoice]], storage[n] - BatteryChargeEfficiency[n] * batterycharge[n] + batterydischarge[n]/BatteryDischargeEfficiency[n] - storageTrials[n, TimeChoice-1, SampleChoice, iter] == 0.0 )
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
    @constraint(model, Pin_Flow_equality[n in HeadNodes[LayerChoice]], p_in[n] - pflow[n-1] == 0.0);

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
    # if TimeChoice == H
    #     @constraint(model, theta == 0 )
    # else
    #     @constraint(model, Cuts[i = 1:length(eCut[LayerChoice, TimeChoice])],
    #     theta  >= eCut[LayerChoice, TimeChoice][i] - sum(ECut[LayerChoice, TimeChoice, n][i]*storage[n] for n in LayerNodes[LayerChoice]) );
    # end
    ### TAKU add
    if TimeChoice == H # no cuts
        @constraint(model, theta == 0.0 )
    else
        if isassigned(eCut[LayerChoice, TimeChoice], OutcomeChoice) # if there exist cuts
            @constraint(model, Cuts[i = 1:length(eCut[LayerChoice, TimeChoice][OutcomeChoice])],
            theta  >= eCut[LayerChoice, TimeChoice][OutcomeChoice][i] - sum(ECut[n,TimeChoice][OutcomeChoice][i] .* storage[n] for n in LayerNodes[LayerChoice]) );
        end
    end
    ### end TAKU add

    ##### Solve the Model
    solve(model)
    println("       Stage Cost value: ", getobjectivevalue(model))
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
        if TimeChoice == H || !isassigned(eCut[LayerChoice, TimeChoice], OutcomeChoice) # TAKU modifided: at the last stage or there is no cut
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
        if TimeChoice == H || !isassigned(eCut[LayerChoice, TimeChoice], OutcomeChoice) #  TAKU modifided: at the last stage or there is no cut
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

    # TAKU add
    # Generate path: THE PATHS NEED TO BE GENERATED BY THE TRUE PROBABILITY (BY USING TransProb)
    SampleScenario = SamplePath(TransProb, K)
    # end TAKU add

    # Solve the First Stage
    println("   ====Forward Pass: solving Layer ",LayerChoice,", Stage ", 1, ", Outcome ",1, ", Iteration ", iter)
    FirstNLDS = NLDS(LayerChoice, 1, 1, 1, iter);
    println("\n")
    # Store The Lower Bound
    LowerBound[LayerChoice, iter] = FirstNLDS.OptimalValue
    # Store the results
    for n in LayerLines[LayerChoice]
        pflowTrials[n, 1, :, iter] = FirstNLDS.pflow[n]*ones(K);
    end
    for n in LayerNodes[LayerChoice]
        storageTrials[n, 1, :, iter] = FirstNLDS.storage[n]*ones(K);
        batterychargeTrials[n, 1, :, iter] = FirstNLDS.batterycharge[n]*ones(K);
        batterydischargeTrials[n, 1, :, iter] = FirstNLDS.batterydischarge[n]*ones(K);
        loadsheddingTrials[n, 1, : , iter] = FirstNLDS.loadshedding[n]*ones(K);
        productionsheddingTrials[n, 1, :, iter] = FirstNLDS.productionshedding[n]*ones(K);
    end
    if LayerChoice == 1
        pgenerationTrials[1:NGenerators, 1, :, iter] = repmat(FirstNLDS.pgeneration,1,K)
    end

    # Solve the remaining stages
    for SampleChoice = 1:K
        for TimeChoice = 2:H
            # Solve the Stage
            # path[TimeChoice,SampleChoice,iter] = rand(1:NLattice[TimeChoice]);  #TAKU: This should be generated by using MC sampling (i.e. Trans prob)
            # TAKU add
            # choose the current scenario
            ScenarioChoice = SampleScenario[LayerChoice,TimeChoice,SampleChoice]
            println("   ====Forward Pass: solving Layer ",LayerChoice,", Stage ", TimeChoice, ", Outcome ",ScenarioChoice, ", MC Sample ",SampleChoice, ", Iteration ",iter )
            # NestedLDS = NLDS(LayerChoice, SampleChoice, TimeChoice, path[TimeChoice, SampleChoice, iter], iter);
            NestedLDS = NLDS(LayerChoice, SampleChoice, TimeChoice, ScenarioChoice, iter);
            println("\n")
            # end TAKU add
            # Store the results
            for n in LayerLines[LayerChoice]
                pflowTrials[n, TimeChoice, SampleChoice, iter] = NestedLDS.pflow[n];
            end
            for n in LayerNodes[LayerChoice]
                batterychargeTrials[n, TimeChoice, SampleChoice, iter] = NestedLDS.batterycharge[n];
                storageTrials[n, TimeChoice, SampleChoice, iter] = NestedLDS.storage[n] ;
                batterydischargeTrials[n, TimeChoice, SampleChoice, iter] = NestedLDS.batterydischarge[n];
                loadsheddingTrials[n, TimeChoice, SampleChoice, iter] = NestedLDS.loadshedding[n] ;
                productionsheddingTrials[n, TimeChoice, SampleChoice, iter] = NestedLDS.productionshedding[n];
            end
            if LayerChoice == 1
                pgenerationTrials[1:NGenerators, TimeChoice, SampleChoice, iter] = NestedLDS.pgeneration ;
            end
        end
    end
end

# Function for the BackwardPass
# function BackwardPass(K, LayerChoice, iter)
#     " K -- Refers to the number of Monte Carlo Samples"
#     " LayerChoice -- Refers to the current Layer"
#     " iter -- Refers to the current iteration"
#     for t = 1:H-1
#         for SampleChoice = 1:K
#             etemp = 0.0;
#             Etemp = Dict{Int64,Float64}(n => 0.0 for n in LayerNodes[LayerChoice]);
#             for OutcomeChoice = 1:NLattice[H+1-t]
#                 # The NLDS is solved
#                 println("   ****Backward Pass: solving layer ",LayerChoice," , Stage ", H-t+1, " ,Outcome ",OutcomeChoice, " ,Sample ",SampleChoice )
#                 NesLDS = NLDS(LayerChoice, SampleChoice, H-t+1, OutcomeChoice, iter);
#                 println("   ****")
#                 # The E coefficients are calculated
#                 # The e coeficients are calculated
#                 for n in HeadNodes[LayerChoice]
#                     etemp = etemp + NesLDS.Pin_fix[n]*p_in_data[H-t+1,OutcomeChoice] ;
#                 end
#                 for n in LayerNodes[LayerChoice]
#                     etemp = etemp + NesLDS.Balance[n]*PNetDemand[n,H-t+1][OutcomeChoice] + NesLDS.StorageMax[n]*BatteryCapacity[n] + NesLDS.BatteryChargeMax[n]*BatteryChargeRate[n] + NesLDS.BatteryDischargeMax[n]*BatteryChargeRate[n];
#                     Etemp[n] = Etemp[n] - TransProb[LayerChoice,H-t][OutcomeChoice]*NesLDS.BatteryDynamics[n];
#                 end
#                 for n in LayerLines[LayerChoice]
#                     etemp = etemp + NesLDS.FlowMax[n]*SLimit[n] + 2*NesLDS.FlowMin[n]*SLimit[n] ;
#                 end
#                 for n in LeafNodes[LayerChoice], m in LeafChildren[LayerChoice,n]
#                     etemp = etemp + NesLDS.Pout_fix[[n, m]]*p_out_data[H-t+1,OutcomeChoice];
#                 end
#                 if LayerChoice == 1
#                     etemp = etemp + sum( NesLDS.GenerationMax[g]*PGenerationMax[g, H-t+1][OutcomeChoice] for g = 1:NGenerators ) - sum( NesLDS.GenerationMin[g]*PGenerationMin[g,H-t+1][OutcomeChoice] for g = 1:NGenerators )
#                 end
#                 if t > 1
#                     etemp = etemp + sum( NesLDS.Cuts[:].*eCut[LayerChoice, H+1-t][:]) ;
#                 end
#                 etemp = TransProb[LayerChoice,H-t][OutcomeChoice]*etemp;
#             end
#             # Here we check if the cut is useful
#             if UsefulCut(LayerChoice, H-t, etemp, Etemp) == 0
#                 for n in LayerNodes[LayerChoice]
#                     push!(ECut[LayerChoice, H-t, n], Etemp[n]);
#                 end
#                 push!(eCut[LayerChoice, H-t], etemp);
#             end
#         end
#     end
# end

# Function for the BackwardPass
function BackwardPass(K, LayerChoice, iter)
    " written by Taku"
    " K -- Refers to the number of Monte Carlo Samples"
    " LayerChoice -- Refers to the current Layer"
    " iter -- Refers to the current iteration"
    for t = 0:H-2 # for loop of OptimalValue
        TimeChoice = H - t;
        for SampleChoice = 1:K  # for loop of Monte Carlo trials
            # store the solutions of NLDS(t,k) for all k (k = OutcomeChoice)
            NesLDS_tk = []

            # for loop for NLDS(t,k) -> solve and store the solution to `NesLDS_tk`
            for OutcomeChoice = 1:NLattice[TimeChoice]
                # The NLDS is solved
                println("   ****Backward Pass: solving Layer ",LayerChoice,", Stage ", TimeChoice, ", Outcome ",OutcomeChoice, ", MC Sample ",SampleChoice, ", Iteration ",iter )
                NesLDS = NLDS(LayerChoice, SampleChoice, TimeChoice, OutcomeChoice, iter);
                println("\n")
                push!(NesLDS_tk, NesLDS);
            end
            # end for loop for NLDS(t,k)

            # for loop for compute the cuts of NLDS(t-1,j) for all j (j = OutcomeChoice_1)
            for OutcomeChoice_1 = 1:NLattice[TimeChoice-1]
                # compute Cut Coefficient: E
                EE = Dict{Int64,Float64}(n => 0.0 for n in LayerNodes[LayerChoice]);
                for n in LayerNodes[LayerChoice]
                    EE[n] = sum(
                        TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k]
                        .*NesLDS_tk[k].BatteryDynamics[n].*(-1) for k = 1:NLattice[TimeChoice]
                    )
                end # end computing Cut Coefficient: E

                # compute Cut right-hand side: e
                ee = 0.0
                ee = sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                        sum(NesLDS_tk[k].BatteryDynamics[n] .* 0.0 for n in LayerNodes[LayerChoice])  # one can ignore b/c this term is always zero
                        +sum(NesLDS_tk[k].Balance[n] .* PNetDemand[n,TimeChoice][k] for n in LayerNodes[LayerChoice])
                        +sum(NesLDS_tk[k].FlowMax[i] .* SLimit[i] for i in LayerLines[LayerChoice])
                        +sum(NesLDS_tk[k].FlowMin[i] .* (SLimit[i] - (-SLimit[i])) for i in LayerLines[LayerChoice])
                        +sum(NesLDS_tk[k].StorageMax[n] .* BatteryCapacity[n] for n in LayerNodes[LayerChoice])
                        +sum(NesLDS_tk[k].BatteryChargeMax[n] .* BatteryChargeRate[n] for n in LayerNodes[LayerChoice])
                        +sum(NesLDS_tk[k].BatteryDischargeMax[n] .* BatteryChargeRate[n] for n in LayerNodes[LayerChoice]))
                        for k = 1:NLattice[TimeChoice] )
                if !isempty(HeadNodes[LayerChoice]) # CAUTION needs to be generalized for multi-layer (p_in_data)
                    ee += sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                        sum(NesLDS_tk[k].Pin_fix[n] .* p_in_data[TimeChoice,k] for n in HeadNodes[LayerChoice])
                        +sum(NesLDS_tk[k].Pin_Flow_equality[n] .* 0.0 for n in HeadNodes[LayerChoice])
                    )
                    for k = 1:NLattice[TimeChoice] ) # one can ignore b/c this term is always zero
                end
                if !isempty(LeafNodes[LayerChoice]) # CAUTION needs to be generalized for multi-layer (p_out_data)
                    for n in LeafNodes[LayerChoice], m in LeafChildren[LayerChoice,n]
                        ee += sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                            sum(NesLDS_tk[k].Pout_fix[[n,m]] .* p_out_data[TimeChoice,k])
                        )
                        for k = 1:NLattice[TimeChoice])
                    end
                end
                if LayerChoice == 1  # for the root node (with generators)
                    ee += sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                    sum(NesLDS_tk[k].GenerationMax[g] .* PGenerationMax[g,TimeChoice][k] for g = 1:NGenerators)
                    +sum(NesLDS_tk[k].GenerationMin[g] .* PGenerationMin[g,TimeChoice][k] for g = 1:NGenerators)
                    )
                    for k = 1:NLattice[TimeChoice] )
                end
                if TimeChoice < H  # if there are cuts
                    ee += sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                        sum(NesLDS_tk[k].Cuts .* eCut[LayerChoice, TimeChoice][k])
                    )
                    for k = 1:NLattice[TimeChoice] )
                end
                # end computing Cut RHS: e

                # check the cut is useful or not
                # if the cut is useful we add it otherwise not
                if isassigned(eCut[LayerChoice,TimeChoice-1], OutcomeChoice_1)  # if there already exist some cuts
                    if UsefulCut(LayerChoice, TimeChoice-1, OutcomeChoice_1, ee, EE) # if the cut is useful -> add it
                        for n in LayerNodes[LayerChoice]
                            push!(ECut[n,TimeChoice-1][OutcomeChoice_1],EE[n])
                        end
                        push!(eCut[LayerChoice,TimeChoice-1][OutcomeChoice_1],ee)
                    end
                else  # if there was no cut -> simply add the cut (maybe at the first iteration)
                    for n in LayerNodes[LayerChoice]
                        ECut[n,TimeChoice-1][OutcomeChoice_1] = [EE[n]]
                    end
                    eCut[LayerChoice,TimeChoice-1][OutcomeChoice_1] = [ee]
                end # end checking/adding the cut
            end # end for loop for compute the cuts of NLDS(t-1,j) for all j
        end # end for loop of Monte Carlo trials
    end # end for loop of time = H ,..., 2
end
