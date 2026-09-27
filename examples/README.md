# Experiments

These scripts replace the former notebook with reproducible, CPU-backed
experiment pipelines. Each script runs a simulation through the `MonteCarloKFP`
package API, prints summary statistics, and saves a PNG under
`examples/artifacts/`.

## Environment setup

In Julia, a directory containing `Project.toml` is an **environment** even when
it is not an installable package and has no `src/` directory.
`examples/Project.toml` exists only to isolate the optional `Plots` dependency
from the numerical package.

From the repository root, instantiate this environment once:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## Running the pipelines

Run any experiment, preferably with multiple Julia threads:

```bash
julia --threads=auto --project=examples examples/annulus_harmonic_measure.jl
julia --threads=auto --project=examples examples/square_dirichlet.jl
julia --threads=auto --project=examples examples/double_well_harmonic_measure.jl
```

The defaults are sized for an interactive run rather than a publication-quality
study. Use the environment variables documented at the top of each script to
increase the trajectory count, grid resolution, or time horizon. For a quick
smoke run, for example:

```bash
MCKFP_TRAJECTORIES=100 MCKFP_MAX_STEPS=100 \
  julia --project=examples examples/annulus_harmonic_measure.jl
```

The `MCKFP_*` variables configure a simulation; they do not install Julia
dependencies.

Kinetic states use package order `(v, x)`. Plots use the conventional axes `x`
horizontally and `v` vertically, so the scripts reorder coordinates only at
presentation time.

The square experiment prescribes data only on the Wiener-regular velocity
boundaries `v = ±1` and the transport-outgoing sides (`x = -1, v < 0` and
`x = 1, v > 0`). Values on the influx and grazing sides (`x = -1, v ≥ 0` and
`x = 1, v ≤ 0`) are one-sided traces estimated by simulation; they are not
assigned Dirichlet values.
