# AntibioticsCellFate-paper

Code for the figures in "Slow trajectories generate divergent cell fates
following antibiotic stress" (Levien, Carvalho, Huynh, Schultz).

## Running

Tested with Julia 1.11.

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/make_all.jl
```

Plots are made with PythonPlot, which installs its own matplotlib the first
time it loads. PDFs go to `figures/`. Each script also runs on its own.

## Contents

`scripts/figN.jl` writes `figures/figN.pdf`, the complete Fig. N of the paper,
schematic panels included. Fonts, colours and figure widths are set in
`src/common.jl`.

- `src/model.jl`: the full model (`Params`, `simulate`).
- `src/common.jl`: packages, paths and plot style shared by the scripts.
- `src/reduced.jl`: nullclines, fixed points and eigenvalues of the reduced
  models.
- `data/schultz_data.csv`: single-cell traces, one row per lineage and frame.
  `x`, `y`, `z` are the TetA, TetR and constitutive reporters, `g` is the
  growth rate and `arr` marks arrested lineages.
