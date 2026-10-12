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

| Script | Output | Paper figure |
| --- | --- | --- |
| `scripts/fig1.jl` | `fig1a.pdf` | 1A |
| `scripts/fig2.jl` | `fig2.pdf` | 2 |
| `scripts/fig3.jl` | `fig3.pdf` | 3 |
| `scripts/fig4.jl` | `fig4_ab.pdf`, `fig4_c.pdf` | 4 |
| `scripts/fig5.jl` | `fig5.pdf` | 5 |
| `scripts/fig6.jl` | `fig6.pdf` | 6 |

Fig. 1B-D are schematics, and Figs. 1 and 4 were assembled from their panels
by hand.

- `src/model.jl`: the full model (`Params`, `simulate`).
- `src/reduced.jl`: nullclines, fixed points and eigenvalues of the reduced
  models.
- `data/schultz_data.csv`: single-cell traces, one row per lineage and frame.
  `x`, `y`, `z` are the TetA, TetR and constitutive reporters, `g` is the
  growth rate and `arr` marks arrested lineages.
