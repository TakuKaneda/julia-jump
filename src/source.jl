# Source file
# contains some functions that are used in policies
# - read csv file of node data and convert it to a proper dataframe
# - convert json file (i.e. stoch params) to array
# - sampling from Trans prob

import StatsBase, JSON

function Read_nodes_csv(path_to_csv)
    "
    read the csv file of data of nodes
    and convert children of nodes to Array
    "
    nodes_df = CSV.read(path_to_csv)

    # convert Children to Array
    s = []
    for i = 1:nrow(nodes_df)
        if typeof(nodes_df[i,:Children]) == Missing
            push!(s,[])
        else
            push!(s,map(parse,split(nodes_df[i,:Children])))
        end
    end
    nodes_df[:Children] = s  # assign
    return nodes_df
end

function ConvertLayerData2Array(path_to_json_LNodes, path_to_json__LLines)
    "
    convert .json LayerNodes/LayerLines data into Array
    Output: -LayerNodes[l] = Array{Int64,1} : ID of nodes in layer l
            -LayerLines[l] = Array{Int64,1} : ID of lines in layer l
    "
    lnodes = String(read(path_to_json_LNodes))
    llines = String(read(path_to_json__LLines))
    lnjson = JSON.parse(lnodes)
    lljson = JSON.parse(llines)

    LayerNodes = Array{Array}(length(lnjson))
    LayerLines = Array{Array}(length(lljson))
    for l=1:length(lnjson)
        LayerNodes[l] = convert(Array{Int64,1},lnjson[string(l)])
        LayerLines[l] = convert(Array{Int64,1},lljson[string(l)])
    end
    return LayerNodes, LayerLines
end

function ConvertHeadLeafNodes2Array(path_to_json_HNodes,path_to_json_LNodes)
    "
    convert .json HeadNodes/LeafNodes data into Array
    Output: -HeadNodes[l] = Array{Int64,1} : IDs of nodes in layer l
            -LeafNodes[l] = Array{Int64,1} : IDs of lines in layer l
            note: elements can be empty array
    "
    hnodes = String(read(path_to_json_HNodes))
    lnodes = String(read(path_to_json_LNodes))
    # lchildren = String(read("../data/two_LeafChildren.json"))
    hnodes_j = JSON.parse(hnodes)
    lnodes_j = JSON.parse(lnodes)
    # lchild_j = JSON.parse(lchildren)

    HeadNodes = Array{Array}(length(hnodes_j))
    LeafNodes = Array{Array}(length(lnodes_j))
    for l=1:length(hnodes_j)
        HeadNodes[l] = convert(Array{Int64,1},hnodes_j[string(l)])
        LeafNodes[l] = convert(Array{Int64,1},lnodes_j[string(l)])
    end
    return HeadNodes, LeafNodes
end

function ConvertLeafChildren2Array(path_to_json, LeafNodes, NNodes)
    "
    convert .json HeadNodes/LeafNodes data into Array
    Output: - LeafChildren[l,n] = Array{Int64,1} : IDs of nodes in layer l
    "
    lchildren = String(read(path_to_json))
    lchild_j = JSON.parse(lchildren)
    NLayers = size(LeafNodes,1)
    LeafChildren = Array{Array}(NLayers,NNodes)
    for l=1:NLayers
        if haskey(lchild_j, string(l))
            dic = lchild_j[string(l)]
            for n in LeafNodes[l]
                if haskey(dic, string(n))
                    LeafChildren[l,n] = convert(Array{Int64,1},dic[string(n)])
                end
            end
        end
    end

    # fill #undef by Int64[]
    for l=1:NLayers
        for n = 1:NNodes
            if ~isassigned(LeafChildren,l,n)
                LeafChildren[l,n] = Int64[]
            end
        end
    end
    return LeafChildren
end


function ConvertPNetDemand2Array(path_to_json)
    "
    convert .json PNetDemand data into Array
    PNetDemand[n,t] = Array{Float64,1} : possible ND value for each outcome
                                         at node n & stage t
       e.g. PNetDemand[2,3][5]: PNetDemand at node 2 (1 in the paper) at outcome 5 of stage 3
    "
    ndstr = String(read(path_to_json)) # read .json file as String
    ndjson = JSON.parse(ndstr)  # convert to Array{Any}

    # convert to Float64
    PNetDemand = Array{Array}(length(ndjson),length(ndjson["1"]))
    for n=1:length(ndjson)
        for t=1:length(ndjson["1"])
            PNetDemand[n,t] = convert(Array{Float64,1},ndjson[string(n-1)][string(t)])
        end
    end
    return PNetDemand
end

function ConvertPGenerationCapacity2Array(path_to_json)
    "
    convert .json PGeneration Capacity data (Max or Min) into Array
    PGenerationMax/Min[g,t] = Array{Float64,1} : possible Capacity value for each outcome
                                         at generator g & stage t
       e.g. PGenerationMax/Min[1,2][3]: PGenerationMax/Min at generator 1 at outcome 3 of stage 2

    "
    capstr = String(read(path_to_json)) # read .json file as String
    json = JSON.parse(capstr);  # convert to Array{Any}
    # convert to Float64
    PGenerationCap = Array{Array}(length(json),length(json["1"]))
    for g=1:length(json)
        for t=1:length(json["1"])
            PGenerationCap[g,t] = convert(Array{Float64,1},json[string(g)][string(t)])
        end
    end
    return PGenerationCap
end


function ConvertTransProb2Array(path_to_json)
    "
    convert TransProb data into Array
    TransProb[l,t] = Array{Float64,2}: Transition probability matrix from
                                       outcome k of stage t to outcome j of
                                       stage t+1 in layer l
      e.g. TransProb[1,2][3,4]: TransProb from outcome 3 of stage 2 to outcome 4 of stage 3 in layer 1
    "
    tpstr = String(read(path_to_json)) # read .json
    tpjson = JSON.parse(tpstr)  # convert to Array{Any}

    # convert to Float64
    TransProb = Array{Array}(length(tpjson),length(tpjson["1"]))
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
    return TransProb
end


function SamplePath(TransProb, NSamples=1)
    "
    Generate sample paths using TransProb.
    return the resulting sample path of outcome in each layer at every stage
    sample_path[l,t,i]: outcome of layer l at stage t of sample i
    "
    sample_path = ones(Int64,size(TransProb)[1],size(TransProb)[2]+1, NSamples)
    for i = 1:NSamples
        for l=1:size(TransProb)[1]
            for t=1:size(TransProb)[2]
                sample_path[l,t+1,i] = StatsBase.sample(1:length(TransProb[l,t][sample_path[l,t,i],:]),
                                    StatsBase.Weights(TransProb[l,t][sample_path[l,t,i],:]))
            end
        end
    end
    return sample_path
end
