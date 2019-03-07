# julia-jump

## Implement Multi-Stage Stochastic Linear Programming

Using Julia `v1.0` with JuMP `v0.19.0`  
**Translation has been done only for SDDP files.**

**Note**: the codes assume you use `Ipopt` (*other solvers, e.g. `Gurobi` had an error*)  

<!-- ### `train-*.jl`: proto-codes for implementing a policy
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
* `*.csv`: contains (deterministic) problem parameters -->

### The stochastic dual dynamic programming (SDDP)
<!-- One can implement the SDDP with running `RunSDDP.jl`. -->
`LoadDataSDDP.jl` loads data for SDDP implementation.
* Usual SDDP: NLDS is built at every time
  - `RunSDDP.jl`: runs the SDDP
  - `SDDP.jl`: functions used for SDDP implementation
  <!-- - `LoadDataSDDP.jl`: loads data for SDDP implementation -->


* Pre-defined SDDP: Build all NLDSs at first
  - `RunSDDPV2.jl`: runs the SDDP
  - `SDDPV2.jl`: functions used for SDDP implementation
  <!-- - `LoadDataSDDP.jl`: loads data for SDDP implementation -->


* Parallel SDDP: Implement SDDP in parallel
  - `RunSDDPParallel.jl`: runs the SDDP
  - `SDDPParallel.jl`: functions used for SDDP implementation
