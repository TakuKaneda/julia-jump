# julia-jump
---

## Translate AMPL to JuMP
**Note: assume you use Gurobi**

### `train-*.jl`: proto-codes for implementing a policy
* `train-perfectforesight.jl`: implements perfect foresight policy for a deterministic problem
* `train-cempc.jl`: implements certainty-equivalent MPC
* `train-sbrmpc.jl`: implements scenario-based robust MPC

### `data/`: only two-layer problem for now
* `*.json`: contain the corresponding stochastic parameters
* `*.csv`: contain (deterministic) problem parameters
