import JSON

## test json
s = "{
    \"1\": [[0.4, 0.6]],
    \"2\": [[0.2, 0.8],[0.5, 0.5]],
    \"3\": [[0.2, 0.8],[0.5, 0.5]]
}"
tp = JSON.parse(s)
## read .json file
ndstr = String(read("data/two_ND.json"))
tpstr = String(read("data/two_TP.json"))
NetDemand = JSON.parse(ndstr)
TransProb = JSON.parse(tpstr)
## sample from TransProb
