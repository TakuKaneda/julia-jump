import JSON

##
lnodes = String(read("../data/two_LayerNodes.json"))
llines = String(read("../data/two_LayerLines.json"))
lnjson = JSON.parse(lnodes)
lljson = JSON.parse(llines)
##
LayerNodes = Array{Array}(length(lnjson))
LayerLines = Array{Array}(length(lljson))
## Convert to array with Int64
for l=1:length(lnjson)
    LayerNodes[l] = convert(Array{Int64,1},lnjson[string(l)])
    LayerLines[l] = convert(Array{Int64,1},lljson[string(l)])
end
## Head/Leaf/Leafchildren nodes
hnodes = String(read("../data/two_HeadNodes.json"))
lnodes = String(read("../data/two_LeafNodes.json"))
lchildren = String(read("../data/two_LeafChildren.json"))
hnodes_j = JSON.parse(hnodes)
lnodes_j = JSON.parse(lnodes)
lchild_j = JSON.parse(lchildren)

##
HeadNodes = Array{Array}(length(hnodes_j))
LeafNodes = Array{Array}(length(lnodes_j))

for l=1:length(hnodes_j)
    HeadNodes[l] = convert(Array{Int64,1},hnodes_j[string(l)])
    LeafNodes[l] = convert(Array{Int64,1},lnodes_j[string(l)])
end
##
LeafChildren = Array{Array}(length(hnodes_j),NNodes)

for l=1:length(hnodes_j)
    if haskey(lchild_j, string(l))
        dic = lchild_j[string(l)]
        for n in LeafNodes[l]
            if haskey(dic, string(n))
                LeafChildren[l,n] = convert(Array{Int64,1},dic[string(n)])
            end
        end
    end
end

for l=1:length(hnodes_j)
    for n = 1:NNodes
        if ~isassigned(LeafChildren,l,n)
            LeafChildren[l,n] = Int64[]
        end
    end
end

# sparse(A)
## test json

s = "{
    \"1\": [[0.4, 0.6]],
    \"2\": [[0.2, 0.8],[0.5, 0.5]],
    \"3\": [[0.2, 0.8],[0.5, 0.5]]
}"
tp = JSON.parse(s)
## read .json file
ndstr = String(read("../data/two_ND.json"))
tpstr = String(read("../data/two_TP.json"))
ndjson = JSON.parse(ndstr)
tpjson = JSON.parse(tpstr)
## convert data into Array
# NetDemand[n,t] = Array{Float64,1} : possible ND value for each outcome
#                                     at node n & stage t
#   e.g. NetDemand[2,3][5]: Netdemand at node 2 (1 in the paper) at outcome 5 of stage 3
# TransProb[l,t] = Array{Float64,2}: Transition probability matrix from
#                                    outcome k of stage t to outcome j of
#                                    stage t+1 in layer l
#   e.g. TransProb[1,2][3,4]: TransProb from outcome 3 of stage 2 to outcome 4 of stage 3 in layer 1
NetDemand = Array{Array}(length(ndjson),length(ndjson["1"]))
TransProb = Array{Array}(length(tpjson),length(tpjson["1"]))
for n=1:length(ndjson)
    for t=1:length(ndjson["1"])
        NetDemand[n,t] = convert(Array{Float64,1},ndjson[string(n-1)][string(t)])
    end
end
for l=1:length(tpjson)
    for t=1:length(tpjson["1"])
        for k=1:length(tpjson[string(l)][string(t+1)])
            if k == 1
                TransProb[l,t] = convert(Array{Float64,1},tpjson[string(l)][string(t+1)][k]).'
            else
                TransProb[l,t] = vcat(TransProb[l,t],convert(Array{Float64,1},tpjson[string(l)][string(t+1)][k]).')
            end
        end
    end
end
## Sampling from TransProb
import StatsBase
sample_path = ones(Int64,size(TransProb)[1],size(TransProb)[2]+1)
for l=1:size(TransProb)[1]
    for t=1:size(TransProb)[2]
        sample_path[l,t+1] = StatsBase.sample(1:length(TransProb[l,t][sample_path[l,t],:]),
                            StatsBase.Weights(TransProb[l,t][sample_path[l,t],:]))
    end
end
print(sample_path)
##
include("test-source.jl")
p = SamplePath(TransProb,2)
## Generator Capacity
pmaxstr = String(read("../data/two_PMax.json"))
pminstr = String(read("../data/two_PMin.json"))
pgmax = JSON.parse(pmaxstr);
pgmin = JSON.parse(pminstr);
##
PGenerationMax = Array{Array}(length(pgmax),length(pgmax["1"]))
PGenerationMin = Array{Array}(length(pgmax),length(pgmax["1"]))
for g=1:length(pgmax)
    for t=1:length(pgmax["1"])
        PGenerationMax[g,t] = convert(Array{Float64,1},pgmax[string(g)][string(t)])
        PGenerationMin[g,t] = convert(Array{Float64,1},pgmin[string(g)][string(t)])
    end
end
