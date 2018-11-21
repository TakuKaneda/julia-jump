include("SDDP.jl")

# Number of Samples
K = 50;

# Number of Iterations
Iterations = 15;

# Definition of the Arrays to store the results
# May not be the best way to store things
pflowTrials = zeros(NLines, H, K, Iterations);
storageTrials = zeros(NNodes, H, K, Iterations);
batterychargeTrials = zeros(NNodes, H, K, Iterations);
batterydischargeTrials = zeros(NNodes, H, K, Iterations);
loadsheddingTrials = zeros(NNodes, H, K, Iterations);
productionsheddingTrials = zeros(NNodes, H, K, Iterations);
pgenerationTrials = zeros(NGenerators, H, K, Iterations);

# To store SDDP result
LowerBound = zeros(NLayers, Iterations) ; # solution of NLDS_{1,1}
SampleCost = zeros(NLayers, K, Iterations); # total cost of each MC sample
MeanCost =  zeros(NLayers, Iterations) ; # mean total cost of MC samples
MeanCostStd = zeros(NLayers, Iterations); # std of mean cost
SDDPTime = zeros(NLayers, 2, Iterations);  # store Forward/BackwardPass time

@time for LayerChoice = 1:NLayers
    for i = 1:Iterations
        tic()
        ForwardPass(K, LayerChoice, i)
        SDDPTime[LayerChoice,1,i] = toc() # time of forward pass
        tic()
        BackwardPass(K, LayerChoice, i)
        SDDPTime[LayerChoice,2,i] = toc() # time of backward pass
    end
end

# Compute the results
for i = 1:Iterations
    for l = 1:NLayers
        for m = 1:K
            if l == 1
                SampleCost[l,m,i] = sum(pgenerationTrials[g, t, m, i]*MargCost[g] for g =1:NGenerators, t = 1:H) + VOLL*sum(loadsheddingTrials[n, t, m, i] for n in LayerNodes[l], t = 1:H)
            else
                SampleCost[l,m,i] = VOLL*sum(loadsheddingTrials[n, t, m, i] for n in LayerNodes[l], t = 1:H);
            end
        end
        MeanCost[l,i] = mean(SampleCost[l,:,i]);
        MeanCostStd[l,i] = std(SampleCost[l,:,i],corrected = false)/sqrt(K)
    end
end

##
# plot the convergence
using Plots
pyplot()
plts = Array{Any}(NLayers)
for l = 1:NLayers
    plts[l] = plot(LowerBound[l,:], title = ("Layer "*string(l)),xaxis = "Iteration",yaxis="Cost (\$)",linecolor = :blue, label="Lower Bound")
    plts[l] = plot!(MeanCost[l,:], linecolor = :red, label="Mean Cost")
    plts[l] = plot!(MeanCost[l,:] + 1.96.*MeanCostStd[l,:], linecolor = :red, linestyle = :dot, label="95% CI")
    plts[l] = plot!(MeanCost[l,:] - 1.96.*MeanCostStd[l,:], linecolor = :red, linestyle = :dot, label="95% CI")
end
plot(plts[1],plts[2],layout=(1,NLayers))
