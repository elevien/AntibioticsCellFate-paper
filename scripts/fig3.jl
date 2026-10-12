# Fig. 3: leading eigenvalue of the 3D model for several γ, and trajectories
# for different initial resources r(0).
include(joinpath(@__DIR__, "..", "src", "common.jl"))

# single gene: u = [x, r, d, λ]
params(; kw...) = Params(; sectors=[:P], A=ones(1, 1), h=3.0, kw...)
p = params()

d_grid = collect(range(1.0, 5000.0, length=4000))
λ_fp   = growth.(d_grid, Ref(p))
D_fp   = D_at_fixed_point.(d_grid, Ref(p))

γs       = [0.0, 0.5, 1.0]
γ_colors = [PythonPlot.cm.viridis(v) for v in range(0.1, 0.9, length=length(γs))]

Ds       = [2000.0, 3000.0, 4000.0]
D_colors = ["tab:blue", "tab:orange", "tab:red"]
r0s      = [0.1, 0.5, 1.0]
r0_ls    = ["-", "--", ":"]

fig, axs = subplots(1, 3, figsize=(10, 4))

# --- leading eigenvalue --------------------------------------------------------
ax = axs[0]
twin = ax.twinx()
twin.plot(D_fp .* UG, λ_fp; color="0.45", lw=1.5, alpha=0.35)
for (γ, col) in zip(γs, γ_colors)
    ax.plot(D_fp .* UG, max_re_3d.(d_grid, γ, Ref(p)); lw=2.0, color=col, label="γ = $γ")
end
ax.axhline(0; color="k", lw=1, ls="--", alpha=0.6)
ax.set_xlabel("External drug \$D\$ (\$\\mu\$g/ml)")
ax.set_ylabel("Max eigenvalue")
ax.set_xlim(0, maximum(D_fp) * UG)
ax.legend(fontsize=11, frameon=false)
twin.tick_params(axis="y", colors="0.3")
twin.set_ylim(0, 1.05)
panel_label(ax, "B"; pad=1.5)

# --- λ(t) and x(t) -------------------------------------------------------------
for (D, col) in zip(Ds, D_colors), (r0, ls) in zip(r0s, r0_ls)
    df = simulate(params(D=D), [0.01, r0, 0.0, p.lambda0]; tmax=30.0)
    lbl = r0 == r0s[1] ? "D = $(D * UG) \$\\mu\$g/ml" : ""
    axs[1].plot(df.time, df.lambda; color=col, lw=1.8, linestyle=ls, label=lbl)
    axs[2].plot(df.time, df.x; color=col, lw=1.8, linestyle=ls)
end
for (r0, ls) in zip(r0s, r0_ls)
    axs[1].plot([], []; color="gray", lw=1.5, linestyle=ls, label="r(0) = $r0")
end
axs[1].set_xlabel("Time")
axs[1].set_ylabel("Growth rate \$\\lambda(t)\$")
axs[1].set_ylim(0, 1.05)
panel_label(axs[1], "C"; pad=1.5)

axs[2].set_xlabel("Time")
axs[2].set_ylabel("Efflux expression \$x(t)\$")
axs[2].set_ylim(0, 1.1)
panel_label(axs[2], "D"; pad=1.5)

fig.tight_layout()
handles, labels = axs[1].get_legend_handles_labels()
fig.legend(handles, labels; fontsize=10, frameon=false, ncol=6, loc="lower center",
           bbox_to_anchor=(0.62, 1.0))
save(fig, "fig3.pdf")
