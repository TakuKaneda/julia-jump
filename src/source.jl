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
    if NSamples == 1
        return sample_path[:,:,1]
    else
        return sample_path
    end
end


function ReadSamplePath(path,H=24)
    "
    Read a sample path from a text file
    Samples[l,t,i]: outcome of layer l at stage t of sample i
    "
    m = readdlm(path)  # read the text file as a matrix
    NLayers = size(m,2)
    NSamples = convert(Int64,size(m,1)/H)
    Samples = zeros(Int64,(NLayers, H, NSamples))
    for i=1:NSamples
        for t = 1:H
            idx = convert(Int64,t+(i-1)*H)
            Samples[:,t,i] = m[idx,:]
        end
    end
    return Samples
end

function ConvertPinPoutData2Array(path_to_p_in,path_to_p_out,NLayers,NNodes,H=24)
    "
    Read the interface flow data
    Output: - p_in_data[l,t][k]: Array{Float64,1}  p_in data at the head node of layer l
                                at stage t and outcome k
            - p_out_data[n,m][t][k]: Array{Float64,1} p_out data from a leaf node n (of layer Node2Layer[n])
                                to the head node m (one of the children) at stage t, outcome k
    "
    # for p_in
    path = String(read(path_to_p_in)) # read .json
    p_in_raw = JSON.parse(path) # raw data

    p_in_data = Array{Array}(NLayers,H) # Array to store data
    [p_in_data[1,t] = Float64[] for t = 1:H]
    for k in keys(p_in_raw)
        for t = 1:H
            p_in_data[parse(Int64,k),t] = convert(Array{Float64,1},p_in_raw[k][string(t)])
        end
    end

    # for p_out[n,m][t][k]
    path = String(read(path_to_p_out)) # read .json
    p_out_raw = JSON.parse(path) # raw data

    p_out_data = Array{Array}(NNodes,NNodes) # Array to store data
    for k1 in keys(p_out_raw)
        n = parse(Int64,k1)
        for k2 in keys(p_out_raw[k1])
            m = parse(Int64,k2)
            p_out_data[n,m] = Array{Any,1}(H)
            for t = 1:H
                p_out_data[n,m][t] = convert(Array{Float64,1},p_out_raw[k1][k2][string(t)])
            end
        end
    end
    return p_in_data, p_out_data
end
