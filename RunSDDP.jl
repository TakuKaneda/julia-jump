# Select the type of problem
problem_size = "multi" # "two" or "multi"

include("SDDP.jl")

# Number of Samples
const K = 5;

# Number of Iterations
const Iterations = 3;

# Definition of the Arrays to store the results
# May not be the best way to store things
pflowTrials = Array{Float64}(NLines, H, K, Iterations);
storageTrials = Array{Float64}(NNodes, H, K, Iterations);
batterychargeTrials = Array{Float64}(NNodes, H, K, Iterations);
batterydischargeTrials = Array{Float64}(NNodes, H, K, Iterations);
loadsheddingTrials = Array{Float64}(NNodes, H, K, Iterations);
productionsheddingTrials = Array{Float64}(NNodes, H, K, Iterations);
pgenerationTrials = Array{Float64}(NGenerators, H, K, Iterations);

# Here we define the arrays to store the cuts
# ECut[:, n, k, t] is a vector of cut-coefficients of node n (`E^l_{n,t,k}` of NLDS of node n at stage t, outcome k)
# eCut[:, k, t, l] is a vector of cut-RHS (`e^l_{t,k}` of NLDS of layer l at stage t, outcome k)
# Note that there are no cuts at stage H
ECut = Array{Float64}(Iterations*K, NNodes, NLattice[2], H-1); # Note: Needs to be generalaized
eCut = Array{Float64}(Iterations*K, NLattice[2], H-1, NLayers); # Note: Needs to be generalaized

# NumCuts[t, l, k] is the number of cuts of layer l at stage t, outcome k
 NumCuts = Array{Int64}(H-1, NLayers, NLattice[2]);

 # We initalize the number of cuts to 0
for i =1:length(NumCuts)
    NumCuts[i] = 0
end

LowerBound = Array{Float64}(NLayers, Iterations) ;
SampleCost = zeros(NLayers, K, Iterations);
MeanCost =  zeros(NLayers, Iterations) ;
MeanCostStd = zeros(NLayers, Iterations); # std of mean cost
SDDPTime = zeros(NLayers, 2, Iterations);  # store Forward/BackwardPass time
#TimeSolve = Array{Float64}(NLayers, K, Iterations);
#TimeOther = Array{Float64}(NLayers, K, Iterations);

@time for LayerChoice = 1:NLayers
    for iter = 1:Iterations
        tic()
        ForwardPass(K, LayerChoice, iter)
        SDDPTime[LayerChoice, 1, iter] = toc()
        # BackwardPass(K, LayerChoice, i)
        tic()
        BackwardPass(K, LayerChoice, iter)
        SDDPTime[LayerChoice, 2, iter] = toc()
    end
end

for iter = 1:Iterations
    for LayerChoice = 1:NLayers
        for SampleChoice = 1:K
            if LayerChoice == 1
                SampleCost[LayerChoice, SampleChoice, iter] = sum(pgenerationTrials[g, TimeChoice, SampleChoice, iter]*MargCost[g] for g =1:NGenerators, TimeChoice = 1:H) + VOLL*sum(loadsheddingTrials[n, TimeChoice, SampleChoice, iter] for n in LayerNodes[LayerChoice], TimeChoice = 1:H);
            else
                SampleCost[LayerChoice, SampleChoice, iter] = VOLL*sum(loadsheddingTrials[n, TimeChoice, SampleChoice, iter] for n in LayerNodes[LayerChoice], TimeChoice = 1:H);
            end
        end
        MeanCost[LayerChoice, iter] = mean(SampleCost[LayerChoice, :, iter]);
        MeanCostStd[LayerChoice, iter] = std(SampleCost[LayerChoice, :, iter],corrected = false)/sqrt(K);
    end
end

##
#plotting the convergence
using Plots
pyplot()
plts = Array{Any}(NLayers)
for l = 1:NLayers
    plts[l] = plot(LowerBound[l,:], title = ("Layer "*string(l)),xaxis = "Iteration",yaxis="Cost (\$)",linecolor = :blue, label="Lower Bound")
    plts[l] = plot!(MeanCost[l,:], linecolor = :red, label="Mean Cost")
    plts[l] = plot!(MeanCost[l,:] + 1.96.*MeanCostStd[l,:], linecolor = :red, linestyle = :dot, label="95% CI")
    plts[l] = plot!(MeanCost[l,:] - 1.96.*MeanCostStd[l,:], linecolor = :red, linestyle = :dot)
end
plot(plts[1],plts[2],layout=(1,NLayers))
