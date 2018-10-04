import StatsBase, JSON


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
