# Fig. 4: trajectories of the three-gene model below and above the critical
# r_y(0) (panels A, B) and for a sweep of r_y(0) (panel C).
include(joinpath(@__DIR__, "..", "src", "common.jl"))

p = Params()
#     x       y    z    r_x    r_y  r_z  d    λ
u0 = [0.0001, 0.1, 1.0, 0.001, 0.6, 0.6, 0.0, 1.0]
with_ry(ry) = [u0[1:4]; ry; u0[6:8]]

# --- A, B ----------------------------------------------------------------------
fig, axs = subplots(figsize=(7, 4), nrows=2, ncols=2, sharex=true)
labels = ["(A.1)" "(B.1)"; "(A.2)" "(B.2)"]
for i in 1:2, j in 1:2
    axs[i-1, j-1].text(0.7, 0.9, labels[i, j], transform=axs[i-1, j-1].transAxes,
        fontsize=15, va="top", ha="left",
        bbox=Dict("facecolor" => "white", "edgecolor" => "none", "alpha" => 0.8,
                  "boxstyle" => "round,pad=0.3"))
end

for (j, ry) in enumerate((0.02, 5.0))
    df = simulate(p, with_ry(ry); tmax=10.0)
    legend = j == 2

    ax = axs[0, j-1]
    ax.axhline(p.Kd / p.D; color="C5", ls="-.", label=L"K_d/D")
    ax.plot(df.time, df.lambda, "k-", label=L"\lambda/\lambda_0")
    ax.plot(df.time, df.d ./ p.D, "k--", label=L"d/D")
    legend && ax.legend(frameon=false, loc=(1.04, 0))

    ax = axs[1, j-1]
    genes = zip(("x", "y", "z"), ("C0", "C3", "C5"))
    for (g, col) in genes
        ax.plot(df.time, df[!, g], col, label=g)
    end
    for (g, col) in genes
        ax.plot(df.time, df[!, "r_"*g], col * "--", label="\$r_$g\$")
    end
    ax.set_xlabel("time (generations)")
    legend && ax.legend(frameon=false, loc=(1.04, 0))
end
axs[0, 0].set_title(L"r_y(0)<r_y^*")
axs[0, 1].set_title(L"r_y(0)>r_y^*")
axs[0, 0].set_ylabel("a.u.")
axs[1, 0].set_ylabel("concentration")
fig.tight_layout()
save(fig, "fig4_ab.pdf")

# --- C -------------------------------------------------------------------------
runs = [simulate(p, with_ry(ry); tmax=8.0) for ry in 0.0001:0.2:4]

fig = figure(figsize=(5, 2))
ax3 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2)
ax2.set_position([0.6, 0.15, 0.4, 0.7])
ax3.view_init(elev=30, azim=20)
ax3.grid(false)
for df in runs
    col = df.lambda[end] < 0.5 ? "C1" : "C0"
    ax3.plot(df.r_y, df.y, df.lambda, "-", color=col, alpha=0.3)
    ax3.plot([df.r_y[end]], [df.y[end]], [df.lambda[end]], ".", color=col)
    ax2.plot(df.y, df.lambda, "-", color=col, alpha=0.3)
    ax2.plot([df.y[end]], [df.lambda[end]], ".", color=col)
end
ax3.set_xlabel(L"$r_y$")
ax3.set_ylabel(L"$y$")
ax3.set_zlabel(L"$\lambda$")
ax2.set_xlabel(L"$y$")
ax2.set_ylabel(L"$\lambda$")
save(fig, "fig4_c.pdf")
