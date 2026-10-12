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

linestyle = Dict(:stable => (ls="-", lw=2.5), :saddle => (ls=":", lw=2.0))

# plot a branch in pieces, one per run of constant stability
function plot_branch(ax, D, y)
    k = 1
    while k <= length(D)
        j = something(findnext(!=(stab[k]), stab, k), length(D) + 1)
        sty = linestyle[stab[k]]
        ax.plot(D[k:j-1], y[k:j-1]; linestyle=sty.ls, linewidth=sty.lw)
        k = j
    end
end

fig, axs = subplots(2, 2, figsize=(8, 6))

# --- A: nullclines -------------------------------------------------------------
ax = axs[0, 0]
D_nc   = [1000.0, 3100.0, 4500.0]
colors = [PythonPlot.cm.plasma(v) for v in (0.2, 0.55, 0.85)]
d_plot = collect(range(0.1, 4999.9, length=3000))

# flow at the intermediate drug level
dv = collect(range(0.0, 5000.0, length=60))
xv = collect(range(0.001, 2.0, length=60))
dg = repeat(dv', length(xv), 1)
xg = repeat(xv, 1, length(dv))
λg = growth.(dg, Ref(p))
dx = λg .* (x_nc.(dg, Ref(p)) .- xg)
dd = p.J .* (D_nc[2] .- dg) .- (p.beta .* xg ./ (p.Kx .+ xg) .+ λg) .* dg
ax.streamplot(dv .* UG, xv, dd .* UG, dx; color="C0", density=1, linewidth=0.5,
              arrowsize=0, zorder=0)

ax.plot(d_plot .* UG, x_nc.(d_plot, Ref(p)); color="k", lw=2.0, label=L"$\dot{x}=0$")
for (D, col) in zip(D_nc, colors)
    ax.plot(d_plot .* UG, d_nc.(d_plot, D, Ref(p)); color=col, lw=1.0, alpha=0.5,
            label="D=$(D * UG) \$\\mu\$g/ml")
    for (d, x, s) in fixed_points_2d(D, d_grid, p)
        if s == :stable
            ax.plot(d * UG, x; marker="o", color=col, ms=8, markeredgewidth=0, zorder=5)
        else
            ax.plot(d * UG, x; marker="o", color=col, ms=8, markerfacecolor="white",
                    markeredgewidth=2, zorder=5)
        end
    end
end
ax.set_xlabel("Intracellular drug \$d\$ (\$\\mu\$g/ml)")
ax.set_ylabel("Efflux expression \$x\$")
ax.set_xlim(10 * UG, 5000 * UG)
ax.set_ylim(0.001, 1.3)
ax.legend(fontsize=11, frameon=false)
panel_label(ax, "A")

# --- B: bifurcation diagram ----------------------------------------------------
ax = axs[0, 1]
plot_branch(ax, D_fp .* UG, λ_fp)
ax.plot([], []; color="gray", ls="-", lw=2.5, label="stable")
ax.plot([], []; color="gray", ls=":", lw=2.0, label="saddle")
ax.set_xlabel("External drug \$D\$ (\$\\mu\$g/ml)")
ax.set_ylabel("Growth rate \$\\lambda^*\$")
ax.set_xlim(0, maximum(D_fp) * UG)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=11, frameon=false)
panel_label(ax, "B")

# --- C: leading eigenvalue (growth rate in grey) -------------------------------
ax = axs[1, 0]
twin = ax.twinx()
plot_branch(ax, D_fp .* UG, maxre)
plot_branch(twin, D_fp .* UG, λ_fp)
for line in twin.get_lines()
    line.set_alpha(0.35); line.set_linewidth(1.5); line.set_color("0.45")
end
ax.axhline(0; color="k", lw=1, ls="--", alpha=0.6)
ax.set_xlabel("External drug \$D\$ (\$\\mu\$g/ml)")
ax.set_ylabel("Max eigenvalue")
ax.set_xlim(0, maximum(D_fp) * UG)
twin.tick_params(axis="y", colors="0.3")
twin.set_ylim(0, 1.05)
panel_label(ax, "C")

# --- D: trajectories -----------------------------------------------------------
ax = axs[1, 1]
for (D, col) in zip([2600.0, 3000.0, 4000.0], ["tab:blue", "tab:orange", "tab:red"])
    t, λ, x = simulate_2d(D, p)
    ax.plot(t, λ; color=col, lw=2.0, linestyle="-", label="D = $(D * UG) \$\\mu\$g/ml")
    ax.plot(t, x; color=col, lw=1.5, linestyle="--")
end
ax.plot([], []; color="k", lw=2.0, linestyle="-", label=L"$\lambda(t)$")
ax.semilogx([], []; color="k", lw=1.5, linestyle="--", label=L"$x(t)$")
ax.set_xlabel("Time")
ax.set_ylabel("Value")
ax.set_ylim(0, 1.1)
ax.legend(fontsize=8, frameon=false, ncol=2, loc="upper right")
panel_label(ax, "D")

fig.tight_layout()
save(fig, "fig2.pdf")
