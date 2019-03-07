using JuMP, Gurobi,  Clp, GLPK, Ipopt, Statistics

# Here we load the SDDP data
include("LoadDataSDDP.jl")

# Definition of the Arrays to store the results
# May not be the best way to store things
pflowTrials = Array{Float64}(undef, NLines, H, K, Iterations);
storageTrials = Array{Float64}(undef, NNodes, H, K, Iterations);
batterychargeTrials = Array{Float64}(undef, NNodes, H, K, Iterations);
batterydischargeTrials = Array{Float64}(undef, NNodes, H, K, Iterations);
loadsheddingTrials = Array{Float64}(undef, NNodes, H, K, Iterations);
productionsheddingTrials = Array{Float64}(undef, NNodes, H, K, Iterations);
pgenerationTrials = Array{Float64}(undef, NGenerators, H, K, Iterations);

# Here we define the arrays to store the cuts
# ECut[:, n, k, t] is a vector of cut-coefficients of node n (`E^l_{n,t,k}` of NLDS of node n at stage t, outcome k)
# eCut[:, k, t, l] is a vector of cut-RHS (`e^l_{t,k}` of NLDS of layer l at stage t, outcome k)
# Note that there are no cuts at stage H
ECut = Array{Float64}(undef, Iterations*K, NNodes, NLattice[2], H-1); # Note: Needs to be generalaized
eCut = Array{Float64}(undef, Iterations*K, NLattice[2], H-1, NLayers); # Note: Needs to be generalaized

# NumCuts[t, l, k] is the number of cuts of layer l at stage t, outcome k
NumCuts = Array{Int64}(undef, H-1, NLayers, NLattice[2]);
# We initalize the number of cuts to 0
for i =1:length(NumCuts)
    NumCuts[i] = 0
end

LowerBound = zeros(NLayers, Iterations) ;
SampleCost = zeros(NLayers, K, Iterations);
MeanCost =  zeros(NLayers, Iterations) ;
MeanCostStd = zeros(NLayers, Iterations); # std of mean cost
SDDPTime = zeros(NLayers, 2, Iterations);  # store Forward/BackwardPass time

TimeSolve = zeros(NLayers, K, Iterations);
TimeOther = zeros(NLayers, K, Iterations);
DefinitionModels = 0.0;

struct MyProblem
  model
  pflow
  storage
  batterycharge
  batterydischarge
  loadshedding
  productionshedding
  p_in
  p_out
  theta
  pgeneration
  previous_storage
  Balance_rootnode
  GenerationMax
  GenerationMin
  Balance
  Balance_headnode
  Balance_leafnode
  Balance_headleafnode
  Pin_fix
  Pout_fix
  Pin_Flow_equality
  FlowMax
  FlowMin
  StorageMax
  BatteryChargeMax
  BatteryDischargeMax
  BatteryDynamics
  Cuts
end

# The struct to Store the solutions
struct Solutions
    " Here we store the decision variables"
    # Saved as Arrays
    pflow#::Array{Float64,1}
    storage#::Array{Float64,1}
    batterycharge#::Array{Float64,1}
    batterydischarge#::Array{Float64,1}
    loadshedding#::Array{Float64,1}
    productionshedding#::Array{Float64,1}
    pgeneration

    " Here we store the dual multipliers"
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
function UsefulCut(LayerChoice, TimeChoice, OutcomeChoice, Etemp, etemp)
    "Written Daniel/Taku"
    "Function that determines wether the obtained cut is repeated
    return  `true` : the cut is useful
            `false`: otherwise
    "
    i = 1;
    boo = true;
    while boo == true && i <= NumCuts[TimeChoice, LayerChoice, OutcomeChoice]
        Difference = sum(abs( ECut[i, LayerNodes[LayerChoice][j], OutcomeChoice, TimeChoice] - Etemp[j]) for j = 1:length(LayerNodes[LayerChoice]) ) # compare E
        Difference += abs(eCut[i, OutcomeChoice, TimeChoice, LayerChoice] - etemp) # compare e
        if Difference < 10e-6
            boo = false; # the cut is repeated i.e. unuseful
        else
            i = i+1; # go to the next cut that will be compared
        end
    end
    return boo
end

function createModel(LayerChoice, TimeChoice, OutcomeChoice)
    # Definition of the model
    # model = Model(with_optimizer(Gurobi.Optimizer, LogToConsole=0)); # ERROR: LoadError: AssertionError: dual <= 0.0 I DONT KNOW WHY
    # model = Model(with_optimizer(Clp.Optimizer, LogLevel=0)); # ERROR: LoadError: AssertionError: dual <= 0.0
    # model = Model(with_optimizer(GLPK.Optimizer, msg_lev=0)); # ERROR: LoadError: AssertionError: dual <= 0.0
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0)); # Ipopt works!

    #model.internalModelLoaded = false
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
    # new variable to be fixed to the previous storage level
    @variable(model, previous_storage[LayerNodes[LayerChoice]]);

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
    # if TimeChoice == 1
    #     @constraint(model,  BatteryDynamics[n in LayerNodes[LayerChoice]], storage[n] - BatteryChargeEfficiency[n] * batterycharge[n] + batterydischarge[n]/BatteryDischargeEfficiency[n] - ini_storage[n] == 0.0 )
    # else
    @constraint(model,  BatteryDynamics[n in LayerNodes[LayerChoice]],
        (storage[n] - BatteryChargeEfficiency[n] * batterycharge[n] + batterydischarge[n]/BatteryDischargeEfficiency[n] == previous_storage[n])
    );
    # end

    # Balance Constraints
    # Balancing - usual nodes
    @constraint(model, Balance[n in UsualNodeSet[LayerChoice] ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m]  for m in Children[n]) - pflow[n-1]  ==  PNetDemand[n,TimeChoice][OutcomeChoice] )  );

    # Balancing - head node
    @constraint(model, Balance_headnode[n in HeadNodeSet[LayerChoice] ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m] for m in Children[n]) + pflow[n-1] ==  PNetDemand[n,TimeChoice][OutcomeChoice] ) );

    # Balancing - leaf node
    @constraint(model, Balance_leafnode[n in LeafNodeSet[LayerChoice] ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m] for m in ChildrenNodesMinusLeafChildren[[LayerChoice,n]]  )
    - pflow[n-1] - sum(p_out[n,j] for j in LeafChildren[LayerChoice,n]) ==  PNetDemand[n,TimeChoice][OutcomeChoice] )  );

    # Balancing - head-leaf node
    @constraint(model, Balance_headleafnode[n in HeadLeafNodeSet[LayerChoice] ],
    (batterydischarge[n] + loadshedding[n] - productionshedding[n] - batterycharge[n]
    + sum(pflow[m] for m in ChildrenNodesMinusLeafChildren[[LayerChoice,n]]  )
    + pflow[n-1] - sum(p_out[n,j] for j in LeafChildren[LayerChoice,n]) ==  PNetDemand[n,TimeChoice][OutcomeChoice] ) );

    # p_in constraint
    @constraint(model, Pin_fix[n in HeadNodes[LayerChoice]], p_in[n] == p_in_data[LayerChoice,TimeChoice][OutcomeChoice]);

    # p_out constraint
    @constraint(model, Pout_fix[n in LeafNodes[LayerChoice], m in LeafChildren[LayerChoice,n]], p_out[n,m] == p_out_data[n,m][TimeChoice][OutcomeChoice]);

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

    if TimeChoice == H # no cuts
        @constraint(model, theta == 0.0 )
    end

    Cuts = ConstraintRef[];

    if LayerChoice == 1
        return MyProblem(model, pflow, storage, batterycharge, batterydischarge, loadshedding, productionshedding, p_in,
        p_out, theta, pgeneration, previous_storage, Balance_rootnode, GenerationMax, GenerationMin, Balance, Balance_headnode,
        Balance_leafnode, Balance_headleafnode, Pin_fix, Pout_fix, Pin_Flow_equality, FlowMax, FlowMin, StorageMax,
        BatteryChargeMax, BatteryDischargeMax, BatteryDynamics, Cuts)

    else
        return MyProblem(model, pflow, storage, batterycharge, batterydischarge, loadshedding, productionshedding, p_in,
        p_out, theta, "no",  previous_storage, "no", "no", "no", Balance, Balance_headnode,
        Balance_leafnode, Balance_headleafnode, Pin_fix, Pout_fix, Pin_Flow_equality, FlowMax, FlowMin, StorageMax,
        BatteryChargeMax, BatteryDischargeMax, BatteryDynamics, Cuts)
    end
end

start = time()
Models = [ [] for LayerChoice = 1:NLayers, TimeChoice = 1:H, OutcomeChoice = 1:NLattice[2] ];
for LayerChoice = 1:NLayers
    for TimeChoice = 1:H
        for OutcomeChoice = 1:NLattice[TimeChoice]
            push!(Models[LayerChoice, TimeChoice, OutcomeChoice], createModel(LayerChoice, TimeChoice, OutcomeChoice))
        end
    end
end
DefinitionModels = time()-start;

# The NLDS algorithm
function NLDS(TimeChoice, SampleChoice, iter, LayerChoice, OutcomeChoice)
    " NLDS(t,k)  t = TimeChoice, k = OutcomeChoice"
    " TimeChoice -- Refers to the time stage t"
    " SampleChoice -- Refers to the index of the sample of the Monte Carlo"
    " iter -- Refers to the current iteration"
    " LayerChoice -- Refers to the current Layer"
    " OutcomeChoice -- Refers to the outcome k of time stage t"
    # Definition of the model
    start = time()
    ModelChoice = Models[LayerChoice, TimeChoice, OutcomeChoice][1];
    # Battery Constraints
    # setRHS has beeen discarded, use 'fix' instead
    if TimeChoice == 1
        for n in LayerNodes[LayerChoice]
            fix(ModelChoice.previous_storage[n],  ini_storage[n]);
        end
    else
        for n in LayerNodes[LayerChoice]
            fix(ModelChoice.previous_storage[n],  storageTrials[n, TimeChoice-1, SampleChoice, iter]);
        end
    end
    TimeOther[LayerChoice, SampleChoice, iter] += time() - start;

    ##### Solve the Model
    start = time()
    optimize!(ModelChoice.model)
    TimeSolve[LayerChoice, SampleChoice, iter] += time() - start;
    #println("       Stage Cost value: ", objective_value(ModelChoice.model))
    ##### Here we return the results
    # Here we translate JuMP.Dict type of dual Multipliers & optimal decisions to Dict Type
    BatteryDynamicsMultipliers =  CreateDictionaryV( push!([], dual.(ModelChoice.BatteryDynamics)) );
    BalanceMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.Balance), dual.(ModelChoice.Balance_headnode), dual.(ModelChoice.Balance_leafnode), dual.(ModelChoice.Balance_headleafnode) ));
    Pin_fixMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.Pin_fix)) );
    Pout_fixMultipliers = CreateDictionaryA( push!([], dual.(ModelChoice.Pout_fix).data) ); # .data is required
    FlowMaxMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.FlowMax)) );
    FlowMinMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.FlowMin)) );
    Pin_Flow_equalityMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.Pin_Flow_equality)) );
    StorageMaxMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.StorageMax)) );
    BatteryChargeMaxMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.BatteryChargeMax)) );
    BatteryDischargeMaxMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.BatteryDischargeMax)) );

    # Return certain data depending on the LayerChoice and TimeChoice
    if LayerChoice == 1
        BalanceMultipliers[1] = dual.(ModelChoice.Balance_rootnode);
        GenerationMaxMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.GenerationMax)) );
        GenerationMinMultipliers = CreateDictionaryV( push!([], dual.(ModelChoice.GenerationMin)) );

        if TimeChoice == H || NumCuts[TimeChoice, LayerChoice, OutcomeChoice] < 1
            return Solutions( value.(ModelChoice.pflow)[:], value.(ModelChoice.storage)[:], value.(ModelChoice.batterycharge)[:],
            value.(ModelChoice.batterydischarge)[:], value.(ModelChoice.loadshedding)[:], value.(ModelChoice.productionshedding)[:],
            value.(ModelChoice.pgeneration)[:], BatteryDynamicsMultipliers, BalanceMultipliers,
            Pin_fixMultipliers, Pout_fixMultipliers, Pin_Flow_equalityMultipliers, FlowMaxMultipliers,
            FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers, BatteryDischargeMaxMultipliers, "NoCuts",
            GenerationMaxMultipliers, GenerationMinMultipliers, objective_value(ModelChoice.model) )
        else
            return Solutions(value.(ModelChoice.pflow)[:], value.(ModelChoice.storage)[:], value.(ModelChoice.batterycharge)[:],
            value.(ModelChoice.batterydischarge)[:], value.(ModelChoice.loadshedding)[:], value.(ModelChoice.productionshedding)[:],
            value.(ModelChoice.pgeneration)[:], BatteryDynamicsMultipliers,
            BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers, Pin_Flow_equalityMultipliers, FlowMaxMultipliers,
            FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers, BatteryDischargeMaxMultipliers,
            dual.(ModelChoice.Cuts), GenerationMaxMultipliers, GenerationMinMultipliers, objective_value(ModelChoice.model) )
        end
    else
        if TimeChoice == H || NumCuts[TimeChoice, LayerChoice, OutcomeChoice] < 1
            return Solutions(value.(ModelChoice.pflow)[:], value.(ModelChoice.storage)[:], value.(ModelChoice.batterycharge)[:],
            value.(ModelChoice.batterydischarge)[:], value.(ModelChoice.loadshedding)[:], value.(ModelChoice.productionshedding)[:],
            "No pgeneration", BatteryDynamicsMultipliers, BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers,
            Pin_Flow_equalityMultipliers, FlowMaxMultipliers, FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers,
            BatteryDischargeMaxMultipliers, "NoCuts", "No pgeneration", "No pgeneration", objective_value(ModelChoice.model) )
        else
            return Solutions(value.(ModelChoice.pflow)[:], value.(ModelChoice.storage)[:], value.(ModelChoice.batterycharge)[:],
            value.(ModelChoice.batterydischarge)[:], value.(ModelChoice.loadshedding)[:], value.(ModelChoice.productionshedding)[:],
            "No pgeneration", BatteryDynamicsMultipliers, BalanceMultipliers, Pin_fixMultipliers, Pout_fixMultipliers,
            Pin_Flow_equalityMultipliers, FlowMaxMultipliers, FlowMinMultipliers, StorageMaxMultipliers, BatteryChargeMaxMultipliers,
            BatteryDischargeMaxMultipliers, dual.(ModelChoice.Cuts), "No pgeneration", "No pgeneration", objective_value(ModelChoice.model) )
        end
    end
end

# This is the inside function for the ForwardPass, helpful to paralelize
function ForwardTrial(iter, LayerChoice, SampleScenario, SampleChoice)
    " LayerChoice -- Refers to the current Layer"
    " iter -- Refers to the current iteration"
    " SampleScenario -- Refers to the array where paths are stored "
    " SampleChoice -- Refers to the current trial solution"
    # Solve the remaining stages
    for TimeChoice = 2:H
        # Solve the Stage
        # choose the current scenario
        OutcomeChoice = SampleScenario[LayerChoice, TimeChoice, SampleChoice]
        #println("   ====Forward Pass: solving layer ",LayerChoice,", Stage ", TimeChoice, ", Outcome ",OutcomeChoice, ", MC Sample ",SampleChoice, ", Iteration ",iter )
        NestedLDS = NLDS(TimeChoice, SampleChoice, iter, LayerChoice, OutcomeChoice);
        #println("\n")
        # Store the results
        pflowTrials[LayerLines[LayerChoice], TimeChoice, SampleChoice, iter] = NestedLDS.pflow;
        batterychargeTrials[LayerNodes[LayerChoice], TimeChoice, SampleChoice, iter] = NestedLDS.batterycharge;
        storageTrials[LayerNodes[LayerChoice], TimeChoice, SampleChoice, iter] = NestedLDS.storage;
        batterydischargeTrials[LayerNodes[LayerChoice], TimeChoice, SampleChoice, iter] = NestedLDS.batterydischarge;
        loadsheddingTrials[LayerNodes[LayerChoice], TimeChoice, SampleChoice, iter] = NestedLDS.loadshedding;
        productionsheddingTrials[LayerNodes[LayerChoice], TimeChoice, SampleChoice, iter] = NestedLDS.productionshedding;
        if LayerChoice == 1
            pgenerationTrials[1:NGenerators, TimeChoice, SampleChoice, iter] = NestedLDS.pgeneration;
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
    #println("   ====Forward Pass: solving layer ",LayerChoice,", Stage ", 1, ", Outcome ",1, ", Iteration ", iter)
    FirstNLDS = NLDS(1, 1, iter, LayerChoice, 1);
    #println("\n")
    # Store The Lower Bound
    LowerBound[LayerChoice, iter] = FirstNLDS.OptimalValue
    # Store the results
    pflowTrials[LayerLines[LayerChoice], 1, :, iter] = repeat(FirstNLDS.pflow,1,K);
    storageTrials[LayerNodes[LayerChoice], 1, :, iter] = repeat(FirstNLDS.storage,1,K);
    batterychargeTrials[LayerNodes[LayerChoice], 1, :, iter] = repeat(FirstNLDS.batterycharge,1,K);
    batterydischargeTrials[LayerNodes[LayerChoice], 1, :, iter] = repeat(FirstNLDS.batterydischarge,1,K);
    loadsheddingTrials[LayerNodes[LayerChoice], 1, :, iter] = repeat(FirstNLDS.loadshedding,1,K);
    productionsheddingTrials[LayerNodes[LayerChoice], 1, :, iter] = repeat(FirstNLDS.productionshedding,1,K);
    if LayerChoice == 1
        pgenerationTrials[1:NGenerators, 1, :, iter] = repeat(FirstNLDS.pgeneration,1,K)
    end
    for SampleChoice = 1:K
        ForwardTrial(iter, LayerChoice, SampleScenario, SampleChoice)
    end
end

# This is the inside function for the BackwardPass, helpful to paralelize
function BackwardTrial(TimeChoice, SampleChoice, iter, LayerChoice, ECutTempStorage, eCutTempStorage)
    " Written by Daniel/Taku"
    " TimeChoice -- Refers to the time stage t"
    " SampleChoice -- Refers to the index of the sample of the Monte Carlo"
    " iter -- Refers to the current iteration"
    " LayerChoice -- Refers to the current Layer"
    " ECutTempStorage -- Array to store temporarily the cut, helpful in the parallel version"
    " eCutTempStorage -- Array to store temporarily the cut, helpful in the parallel version"

    NesLDS_tk = [];
    # for loop for NLDS(t,k) -> solve and store the solution to `NesLDS_tk`
    for OutcomeChoice = 1:NLattice[TimeChoice]
        # The NLDS is solved
        #println("   ****Backward Pass: solving layer ",LayerChoice,", Stage ", TimeChoice, ", Outcome ",OutcomeChoice, ", MC Sample ",SampleChoice, ", Iteration ",iter )
        NesLDS = NLDS(TimeChoice, SampleChoice, iter, LayerChoice, OutcomeChoice)
        #println("\n")
        push!(NesLDS_tk, NesLDS);
    end
    # end for loop for NLDS(t,k)

    # for loop for compute the cuts of NLDS(t-1,j) for all j (j = OutcomeChoice_1)
    for OutcomeChoice_1 = 1:NLattice[TimeChoice-1]
        # compute Cut Coefficient: E
        EE = Array{Float64}(undef,length(LayerNodes[LayerChoice]));
        for i = 1:length(LayerNodes[LayerChoice])
            n = LayerNodes[LayerChoice][i];
            EE[i] = sum(
                TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k]
                .*NesLDS_tk[k].BatteryDynamics[n].*(-1) for k = 1:NLattice[TimeChoice]
            )
        end # end computing Cut Coefficient: E

        # compute Cut right-hand side: e
        ee = 0.0
        ee = sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                +sum(NesLDS_tk[k].Balance[n] .* PNetDemand[n,TimeChoice][k] for n in LayerNodes[LayerChoice])
                +sum(NesLDS_tk[k].FlowMax[i] .* SLimit[i] for i in LayerLines[LayerChoice])
                +sum(NesLDS_tk[k].FlowMin[i] .* (SLimit[i] - (-SLimit[i])) for i in LayerLines[LayerChoice])
                +sum(NesLDS_tk[k].StorageMax[n] .* BatteryCapacity[n] for n in LayerNodes[LayerChoice])
                +sum(NesLDS_tk[k].BatteryChargeMax[n] .* BatteryChargeRate[n] for n in LayerNodes[LayerChoice])
                +sum(NesLDS_tk[k].BatteryDischargeMax[n] .* BatteryChargeRate[n] for n in LayerNodes[LayerChoice]))
                for k = 1:NLattice[TimeChoice] )
        if !isempty(HeadNodes[LayerChoice])
            ee += sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] *
                sum(NesLDS_tk[k].Pin_fix[n] .* p_in_data[LayerChoice,TimeChoice][k] for n in HeadNodes[LayerChoice])
                for k = 1:NLattice[TimeChoice] )
        end
        if !isempty(LeafNodes[LayerChoice])
            for n in LeafNodes[LayerChoice], m in LeafChildren[LayerChoice,n]
                ee += sum(TransProb[LayerChoice,TimeChoice-1][OutcomeChoice_1,k] * (
                    sum(NesLDS_tk[k].Pout_fix[[n,m]] .* p_out_data[n,m][TimeChoice][k])
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
                sum(NesLDS_tk[k].Cuts .* eCut[1:NumCuts[TimeChoice, LayerChoice, k], k, TimeChoice, LayerChoice]
                 ) ) for k = 1:NLattice[TimeChoice] )
        end
        # Here we store the Cuts
        for i = 1:length(LayerNodes[LayerChoice])
            ECutTempStorage[SampleChoice, OutcomeChoice_1, i] = EE[i];
        end
        eCutTempStorage[SampleChoice, OutcomeChoice_1] = ee;
    end # end computing cuts
end

# Function for the BackwardPass
function BackwardPass(K, LayerChoice, iter)
    " K -- Refers to the number of Monte Carlo Samples"
    " LayerChoice -- Refers to the current Layer"
    " iter -- Refers to the current iteration"

    for t = 0:H-2 # for loop of time Stage
        # Define the current time stage
        TimeChoice = H - t;
        # Define the temporal arrays to store the cuts, useful for the parallel version
        ECutTempStorage = Array{Float64}(undef,K, NLattice[TimeChoice], length(LayerNodes[LayerChoice]));
        eCutTempStorage = Array{Float64}(undef,K, NLattice[TimeChoice]);

        for SampleChoice = 1:K # for loop for the trials
            BackwardTrial(TimeChoice, SampleChoice, iter, LayerChoice, ECutTempStorage, eCutTempStorage)
        end

        # Check the cut is useful or not
        # if the cut is useful we add it otherwise not

        for SampleChoice = 1:K
            for OutcomeChoice_1 = 1:NLattice[TimeChoice-1]
                Etemp = ECutTempStorage[SampleChoice, OutcomeChoice_1, :];
                etemp = eCutTempStorage[SampleChoice, OutcomeChoice_1];
                if NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1] >= 1  # if there already exist some cuts
                    if UsefulCut(LayerChoice, TimeChoice-1, OutcomeChoice_1, Etemp, etemp) # if the cut is useful -> add it
                        NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1] += 1;
                        for i = 1:length(LayerNodes[LayerChoice])
                            n = LayerNodes[LayerChoice][i];
                            ECut[NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1], n, OutcomeChoice_1, TimeChoice-1] = Etemp[i];
                        end
                        eCut[NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1], OutcomeChoice_1, TimeChoice-1, LayerChoice] = etemp;

                        start = time();
                        ModelChoice = Models[LayerChoice, TimeChoice-1, OutcomeChoice_1][1];
                        push!(ModelChoice.Cuts, @constraint(ModelChoice.model, ModelChoice.theta >= etemp - sum(ECut[NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1], n, OutcomeChoice_1, TimeChoice-1] * ModelChoice.storage[n] for n in LayerNodes[LayerChoice]) ));
                        TimeOther[LayerChoice, SampleChoice, iter] += time() - start;
                    end
                else  # if there was no cut -> simply add the cut (maybe at the first iteration)
                    NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1] += 1;
                    for i = 1:length(LayerNodes[LayerChoice])
                        n = LayerNodes[LayerChoice][i];
                        ECut[NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1], n, OutcomeChoice_1, TimeChoice-1] = Etemp[i];
                    end
                    eCut[NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1], OutcomeChoice_1, TimeChoice-1, LayerChoice] = etemp;
                    start = time();
                    ModelChoice = Models[LayerChoice, TimeChoice-1, OutcomeChoice_1][1];
                    push!(ModelChoice.Cuts, @constraint(ModelChoice.model, ModelChoice.theta >= etemp - sum(ECut[NumCuts[TimeChoice-1, LayerChoice, OutcomeChoice_1], n, OutcomeChoice_1, TimeChoice-1] * ModelChoice.storage[n] for n in LayerNodes[LayerChoice]) ));
                    TimeOther[LayerChoice, SampleChoice, iter] += time() - start;
                end
            end
        end # end checking/adding the cut
    end # end for loop of time = H ,..., 2
end
