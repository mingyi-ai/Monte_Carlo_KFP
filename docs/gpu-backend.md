# GPU backend

## Current support

MonteCarloKFP currently supports CPU execution only. GPU packages are not
dependencies, and there is no supported accelerator backend or backend-selection
API at this time.

A future GPU implementation should be exposed through the package API rather
than through a hardware-specific notebook. It should remain optional, preserve
the dynamics and domain semantics of `simulate`, use an appropriate reproducible
device RNG, and be validated statistically against the CPU implementation. The
backend interface will be designed when a maintained second implementation
exists.

## Legacy Metal prototype

The initial repository contained a naive Apple Metal prototype in
`[Julia][Metal]Monte-Carlo-KFP.ipynb`, `MCKernel.jl`, and `RandomMetal.jl`.
Experiment orchestration, plotting, device allocation, kernels, and
random-number generation were coupled directly in the notebook workflow. Those
files are no longer part of the working tree, but remain available in the Git
history (beginning with the initial commit `6ca4d5e`).

The experiment portions were replaced by the scripts in
[`examples/`](../examples/), which use the supported CPU package API and save
plots for remote/headless runs. Metal may be supported again in the future
through the same package-facing workflow, but the historical notebook code is
not a supported backend.
