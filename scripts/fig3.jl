# Fig. 3: leading eigenvalue of the 3D model for several γ, and trajectories
# for different initial resources r(0).
include(joinpath(@__DIR__, "..", "src", "common.jl"))

# single gene: u = [x, r, d, λ]
params(; kw...) = Params(; sectors=[:P], A=ones(1, 1), h=3.0, kw...)
p = params()

d_grid = collect(range(1.0, 5000.0, length=4000))
D_fp   = D_at_fixed_point.(d_grid, Ref(p))

γs       = [0.0, 0.5, 1.0]
γ_colors = [PythonPlot.cm.viridis(v) for v in (0.0, 0.45, 0.8)]

Ds    = [2000.0, 3000.0, 4000.0]
r0s   = [0.1, 0.5, 1.0]
r0_ls = ["-", "--", ":"]

fig, axs = subplots(1, 3, figsize=(WIDTH, 3.3))

# --- A: leading eigenvalue -----------------------------------------------------
ax = axs[0]
for (γ, col) in zip(γs, γ_colors)
    ax.plot(D_fp .* UG, max_re_3d.(d_grid, γ, Ref(p)); color=col, label="\$\\gamma=$γ\$")
end
ax.axhline(0; color="0.5", lw=1, ls="--", zorder=0)
ax.set_xlabel(DRUG_LABEL)
ax.set_ylabel("Max eigenvalue")
ax.set_xlim(0, 5)
ax.set_ylim(-1.05, 0.75)
ax.legend(loc="upper left", handlelength=1.0, ncol=2)
panel_label(ax, "A")

# --- B, C: λ(t) and x(t) -------------------------------------------------------
for (D, col) in zip(Ds, DRUG_COLORS), (r0, ls) in zip(r0s, r0_ls)
    df = simulate(params(D=D), [0.01, r0, 0.0, p.lambda0]; tmax=30.0)
    axs[1].plot(df.time, df.lambda; color=col, linestyle=ls)
    axs[2].plot(df.time, df.x; color=col, linestyle=ls)
end
axs[1].set_ylabel("Growth rate \$\\lambda\$")
axs[2].set_ylabel("Efflux expression \$x\$")
for (ax, s) in zip((axs[1], axs[2]), ("B", "C"))
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, 0.5, 1])
    panel_label(ax, s)
end

# colour codes D (top row of the legend), line style codes r(0) (bottom row)
line(; kw...) = PythonPlot.matplotlib.lines.Line2D([], []; kw...)
handles = []
for i in 1:3
    push!(handles, line(color=DRUG_COLORS[i], label="\$D=$(Ds[i] * UG)\$"))
    push!(handles, line(color="k", ls=r0_ls[i], label="\$r(0)=$(r0s[i])\$"))
end
fig.legend(handles=handles, loc="outside upper right", ncol=3, handlelength=1.2)

save(fig, "fig3.pdf")
