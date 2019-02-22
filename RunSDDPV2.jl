# Select the type of problem
problem_size = "multi" # "two" or "multi"

# Number of Samples
const K = 3;

# Number of Iterations
const Iterations = 5;

include("SDDPV2.jl")

for LayerChoice = 1:6
    for iter = 1:Iterations
        tic()
        ForwardPass(K, LayerChoice, iter)
        SDDPTime[LayerChoice, 1, iter] = toq();
        # BackwardPass(K, LayerChoice, i)
        tic()
        BackwardPass(K, LayerChoice, iter)
        SDDPTime[LayerChoice, 2, iter] = toq();
    end
end

for iter = 1:Iterations
    for LayerChoice = 1:6
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
plts = Array{Any}(6)
for l = 1:6
    plts[l] = plot(LowerBound[l,:], title = ("Layer "*string(l)),xaxis = "Iteration",yaxis="Cost (\$)",linecolor = :blue, label="Lower Bound")
    plts[l] = plot!(MeanCost[l,:], linecolor = :red, label="Mean Cost")
    plts[l] = plot!(MeanCost[l,:] + 1.96.*MeanCostStd[l,:], linecolor = :red, linestyle = :dot, label="95% CI")
    plts[l] = plot!(MeanCost[l,:] - 1.96.*MeanCostStd[l,:], linecolor = :red, linestyle = :dot)
end
plot(plts[1],plts[2],plts[3],plts[4],plts[5],plts[6],layout=(2,3))
