# Fig. 2: nullclines, bifurcation diagram, leading eigenvalue and sample
# trajectories of the 2D model.
include(joinpath(@__DIR__, "..", "src", "common.jl"))

p = Params(h=3.0)

# fixed-point branch, parametrised by internal drug level
d_grid = collect(range(1.0, 5000.0, length=4000))
λ_fp   = growth.(d_grid, Ref(p))
D_fp   = D_at_fixed_point.(d_grid, Ref(p))
stab   = stability_2d.(d_grid, Ref(p))
maxre  = max_re_2d.(d_grid, Ref(p))

linestyle = Dict(:stable => "-", :saddle => ":")

# plot a branch in pieces, one per run of constant stability
function plot_branch(ax, D, y; kw...)
    k = 1
    while k <= length(D)
        j = something(findnext(!=(stab[k]), stab, k), length(D) + 1)
        ax.plot(D[k:j-1], y[k:j-1]; linestyle=linestyle[stab[k]], kw...)
        k = j
    end
end

function stability_legend(ax; kw...)
    ax.plot([], []; color="k", ls="-", label="stable")
    ax.plot([], []; color="k", ls=":", label="saddle")
    ax.legend(; kw...)
end

fig, axs = subplots(2, 2, figsize=(WIDTH, 5.4))

# --- A: nullclines -------------------------------------------------------------
ax = axs[0, 0]
D_nc   = [1000.0, 3100.0, 4500.0]
d_plot = collect(range(0.1, 4999.9, length=3000))

# flow at the intermediate drug level
dv = collect(range(0.0, 5000.0, length=60))
xv = collect(range(0.001, 2.0, length=60))
dg = repeat(dv', length(xv), 1)
xg = repeat(xv, 1, length(dv))
λg = growth.(dg, Ref(p))
dx = λg .* (x_nc.(dg, Ref(p)) .- xg)
dd = p.J .* (D_nc[2] .- dg) .- (p.beta .* xg ./ (p.Kx .+ xg) .+ λg) .* dg
ax.streamplot(dv .* UG, xv, dd .* UG, dx; color="0.8", density=0.8, linewidth=0.6,
              arrowsize=0.7, zorder=0)

ax.plot(d_plot .* UG, x_nc.(d_plot, Ref(p)); color="k", lw=2, label="\$\\dot{x}=0\$")
for (D, col) in zip(D_nc, DRUG_COLORS)
    ax.plot(d_plot .* UG, d_nc.(d_plot, D, Ref(p)); color=col,
            label="\$D=$(D * UG)\$")
    for (d, x, s) in fixed_points_2d(D, d_grid, p)
        ax.plot(d * UG, x; marker="o", color=col, ms=7, zorder=5,
                markerfacecolor=(s == :stable ? col : "white"), markeredgewidth=1.5)
    end
end
ax.set_xlabel("Internal drug \$d\$ (µg/ml)")
ax.set_ylabel("Efflux expression \$x\$")
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.3)
ax.legend(loc="upper right", handlelength=1.0)
panel_label(ax, "A")

# --- B: bifurcation diagram ----------------------------------------------------
ax = axs[0, 1]
plot_branch(ax, D_fp .* UG, λ_fp; color="k", lw=2)
stability_legend(ax; loc="lower left")
ax.set_xlabel(DRUG_LABEL)
ax.set_ylabel("Growth rate \$\\lambda^*\$")
ax.set_xlim(0, 5)
ax.set_ylim(0, 1.05)
ax.set_yticks([0, 0.5, 1])
panel_label(ax, "B")

# --- C: leading eigenvalue -----------------------------------------------------
ax = axs[1, 0]
plot_branch(ax, D_fp .* UG, maxre; color="k", lw=2)
ax.axhline(0; color="0.5", lw=1, ls="--", zorder=0)
ax.set_xlabel(DRUG_LABEL)
ax.set_ylabel("Max eigenvalue")
ax.set_xlim(0, 5)
panel_label(ax, "C")

# --- D: trajectories -----------------------------------------------------------
ax = axs[1, 1]
for (D, col) in zip([2600.0, 3000.0, 4000.0], DRUG_COLORS)
    t, λ, x = simulate_2d(D, p)
    ax.semilogx(t, λ; color=col, ls="-", label="\$D=$(D * UG)\$")
    ax.semilogx(t, x; color=col, ls="--")
end
ax.set_xlim(1e-2, 30)
ax.set_ylim(0, 1.05)
ax.set_yticks([0, 0.5, 1])
ax.set_xlabel(TIME_LABEL)
ax.set_ylabel("\$\\lambda\$ (solid), \$x\$ (dashed)")
ax.legend(loc="center left", handlelength=1.0)
panel_label(ax, "D")

save(fig, "fig2.pdf")
