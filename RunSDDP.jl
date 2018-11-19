include("SDDP.jl")

# Number of Samples
K = 50;

# Number of Iterations
Iterations = 15;

# Definition of the Arrays to store the results
# May not be the best way to store things
pflowTrials = zeros(NLayers, NLines, H, K, Iterations);
storageTrials = zeros(NLayers, NNodes, H, K, Iterations);
batterychargeTrials = zeros(NLayers, NNodes, H, K, Iterations);
batterydischargeTrials = zeros(NLayers, NNodes, H, K, Iterations);
loadsheddingTrials = zeros(NLayers, NNodes, H, K, Iterations);
productionsheddingTrials = zeros(NLayers, NNodes, H, K, Iterations);
pgenerationTrials = zeros(NGenerators, H, K, Iterations);
path = Array{Int64,3}(H,K,Iterations);

eCut = [Float64[] for LayerChoice = 1:2, TimeChoice = 1:H-1];
ECut = [Float64[] for LayerChoice = 1:2, TimeChoice = 1:H-1, n=1:NNodes ];

LowerBound = zeros(NLayers, Iterations) ;
MeanCost =  zeros(NLayers, Iterations) ;

@time for LayerChoice = 1:NLayers
    for i = 1:Iterations
        ForwardPass(K, LayerChoice, i)
        BackwardPass(K, LayerChoice, i)
    end
end

for i = 1:Iterations
    for SampleChoice = 1:K
        for TimeChoice = 1:H
            MeanCost[1, i] = MeanCost[1, i] + sum(pgenerationTrials[g, TimeChoice, SampleChoice, i]*MargCost[g] for g =1:3) + VOLL*sum(loadsheddingTrials[1, n, TimeChoice, SampleChoice, i] for n in LayerNodes[1]);
            MeanCost[2, i] = MeanCost[2, i] + VOLL*sum(loadsheddingTrials[2, n, TimeChoice, SampleChoice, i] for n in LayerNodes[2]);
        end
    end
end

MeanCost=MeanCost/K
