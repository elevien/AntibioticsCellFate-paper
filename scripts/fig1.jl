# Fig. 1: experimental trajectories (A) and schematics (B-D).
include(joinpath(@__DIR__, "..", "src", "common.jl"))
Random.seed!(3)

fig = figure(figsize=(WIDTH, 4.6))
axs = fig.subplot_mosaic([["A", "A", "B1", "B2"], ["C", "C", "D", "D"]])

# draw x and y axes as arrows, without ticks
function arrow_axes(ax)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.plot(1, 0, ">k", transform=ax.transAxes, clip_on=false, ms=5)
    ax.plot(0, 1, "^k", transform=ax.transAxes, clip_on=false, ms=5)
end

# --- A: growth rate vs TetR for each lineage ------------------------------------
ax = axs["A"]
data = CSV.read(DATA, DataFrame)
for (arr, col) in ((0, RECOVERED), (1, ARRESTED))
    for lin in unique(data[data.arr .== arr, :lineage])
        d = data[(data.lineage .== lin) .& (data.time .< 50), :]
        ax.plot(d.y, d.g, "-", color=col, alpha=0.3, lw=0.5)
        ax.plot([d.y[1]], [d.g[1]], "+", color=col)
        ax.plot([d.y[end]], [d.g[end]], ".", color=col, ms=8)
    end
end
ax.axhline(0, color="k", ls="--", lw=1)
ax.plot([], [], "-", color=RECOVERED, label="recovered")
ax.plot([], [], "-", color=ARRESTED, label="arrested")
ax.legend(loc="upper right", ncol=2)
ax.set_ylim(-0.7, 2.3)
ax.set_yticks([0, 1])
ax.set_xlabel("TetR (a.u.)")
ax.set_ylabel("Growth rate \$\\lambda/\\langle\\lambda_0\\rangle\$")
panel_label(ax, "A")

# --- B: two routes to heterogeneous fates ---------------------------------------
switch(x) = 2 / (1 + exp(-20 * (x - 1)))
f_fast(x, y) = (y - x, switch(x) - y)
f_slow(x, y) = (y * (y - x), x^0.1 * (switch(x) - y))

function trajectory(f, x, y; T, noise=0.0, dt=0.01)
    xs, ys = [x], [y]
    for _ in 1:round(Int, T / dt)
        dx, dy = f(x, y)
        x = max(x + dt * dx + noise * sqrt(dt) * randn(), 0.0)
        y = max(y + dt * dy + noise * sqrt(dt) * randn(), 0.0)
        push!(xs, x); push!(ys, y)
    end
    return xs, ys
end

function plot_fates(ax, f, starts; kw...)
    for (x0, y0) in starts
        xs, ys = trajectory(f, x0, y0; kw...)
        col = xs[end] > 1 ? RECOVERED : ARRESTED
        ax.plot(xs, ys, color=col, lw=0.8)
        ax.plot([xs[end]], [ys[end]], ".", color=col, ms=8)
    end
    ax.set_xlim(-0.1, 2.6)
    ax.set_ylim(-0.1, 2.6)
    arrow_axes(ax)
end

# fast switch with noise: cells hop between two stable states
plot_fates(axs["B1"], f_fast, [(0.3, 0.3), (1.7, 1.7)]; T=5.0, noise=0.25)
axs["B1"].set_xlabel("Bistability")
# no noise: fates are set by the initial condition, and arrested cells crawl
plot_fates(axs["B2"], f_slow, [(0.0, y0) for y0 in range(0.3, 2.4, length=7)]; T=30.0)
axs["B2"].set_xlabel("Slow time-scale")
panel_label(axs["B1"], "B")

# --- C: model diagram ------------------------------------------------------------
ax = axs["C"]
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
red, green = "tab:red", "tab:green"
arrow(a, b; color="k", style="->") = ax.annotate("", xy=b, xytext=a,
    arrowprops=Dict("arrowstyle" => style, "color" => color, "lw" => 1.2,
                    "shrinkA" => 0, "shrinkB" => 0))

# membrane, with an efflux pump sitting in it
θ = range(-0.75, 0.75, length=50)
ax.plot(0.15 .+ 0.75 .* cos.(θ), 0.55 .+ 0.6 .* sin.(θ), "k", lw=3, solid_capstyle="round")
ax.plot([0.9], [0.55], "o", color=green, ms=9)

# drug molecules inside and outside the cell
ax.plot([0.60, 0.66, 0.72, 0.63, 0.75, 0.82], [0.62, 0.72, 0.60, 0.50, 0.46, 0.64], "^",
        mfc="none", color=red, ms=6)
ax.plot([0.97, 0.99], [0.66, 0.40], "^", mfc="none", color=red, ms=6)
ax.text(0.58, 0.84, "Drug \$d\$", color=red, ha="center")

# drug inhibits growth
ax.text(0.12, 0.60, "\$\\lambda\$", ha="center", va="center",
        bbox=Dict("boxstyle" => "circle,pad=0.3", "fc" => "none", "ec" => "k", "ls" => "--"))
arrow((0.56, 0.60), (0.24, 0.60); color=red, style="|-|,widthA=0,widthB=0.4")
ax.text(0.0, 0.80, "Growth rate", ha="left")

# growth sets expression of the pump
ax.plot([0.03, 0.80], [0.05, 0.05], "k", lw=1)
arrow((0.12, 0.46), (0.12, 0.08))
ax.plot([0.16, 0.16], [0.05, 0.14], "k", lw=1.2)
arrow((0.16, 0.14), (0.28, 0.14))
arrow((0.48, 0.05), (0.48, 0.20); color=green)
ax.plot([0.44, 0.52, 0.47], [0.26, 0.30, 0.36], "o", color=green, ms=7)
ax.text(0.60, 0.26, "Efflux\npump \$x\$", color=green, ha="left", va="center")
panel_label(ax, "C")

# --- D: growth laws ---------------------------------------------------------------
ax = axs["D"]
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
arrow_axes(ax)
arrow(a, b; kw...) = ax.annotate("", xy=b, xytext=a,
    arrowprops=Dict("arrowstyle" => "->", "lw" => 1.2, "shrinkA" => 0, "shrinkB" => 0, kw...))

# nutrient quality moves cells up this line; drug moves them along the others
nutrient(t) = (0.25 + 0.7t, 0.08 + 0.5t)
arrow(nutrient(0), nutrient(1))
ax.text(0.80, 0.34, "+nutrient", rotation=33, ha="center", va="center")
for (t, col) in zip((0.12, 0.35, 0.60, 0.85), ("tab:red", "navy", "tab:green", "k"))
    x0, y0 = nutrient(t)
    x1, y1 = 0.06 + 0.04t, 0.95
    ax.plot([x0, x1], [y0, y1], color=col, lw=1.2)
    s = [0.0, 0.3, 0.6]
    ax.plot(x0 .+ s .* (x1 - x0), y0 .+ s .* (y1 - y0), "o", color=col, ms=5)
end
ax.text(0.68, 0.80, "+drug", rotation=-31, ha="center", va="center")
arrow((0.62, 0.72), (0.38, 0.86))
ax.set_xlabel("Growth rate")
ax.set_ylabel("Ribosome content")
panel_label(ax, "D")

save(fig, "fig1.pdf")
