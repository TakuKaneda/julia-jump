# julia-jump

## Translate AMPL to JuMP
**Note: the codes assume you use Gurobi**

### `train-*.jl`: proto-codes for implementing a policy
* `train-perfectforesight.jl`: implements perfect foresight policy for a deterministic problem
* `train-cempc.jl`: implements certainty-equivalent MPC
* `train-sbrmpc.jl`: implements scenario-based robust MPC

### `data/`: two-layer or multi-layer problem of my master thesis
* `*.json`: contains the corresponding stochastic parameters
* `*.csv`: contains (deterministic) problem parameters
