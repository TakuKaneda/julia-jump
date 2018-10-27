# julia-jump

## Translate AMPL to JuMP
**Note: the codes assume you use Gurobi**

### `train-*.jl`: proto-codes for implementing a policy
* `train-perfectforesight.jl`: implements perfect foresight policy for a deterministic problem
* `train-cempc.jl`: implements certainty-equivalent MPC
* `train-sbrmpc.jl`: implements scenario-based robust MPC

### `parallel-*.jl`: main-codes for implementing a policy in parallel
`prarallel/`: dir which contains the MPC functions or data loading functions
* `parallel-run-perfectforesight.jl`: implements perfect foresight policy for a deterministic problem in parallel
* `parallel-run-cempc.jl`: implements certainty-equivalent MPC in parallel
* `parallel-run-sbrmpc.jl`: implements scenario-based robust MPC in parallel

### `src/`: source codes
* `source.jl`: functions used in the policies such as reading data files

### `data/`: two-layer or multi-layer problem of my master thesis
* `*.json`: contains the corresponding stochastic parameters
* `*.csv`: contains (deterministic) problem parameters
