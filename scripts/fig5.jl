# Fig. 5: experimental lineages (A-C) against simulations over a range of
# initial TetR resources r_y(0) (D-F).
include(joinpath(@__DIR__, "..", "src", "common.jl"))
Random.seed!(1234)

# per-lineage final growth rate and initial velocity of y
summarise(df) = combine(groupby(df, :lineage),
    :lambda => (v -> mean(v[end-10:end])) => :lambda,
    :y => (v -> mean(diff(v[1:4]))) => :yv)

# --- experiment ----------------------------------------------------------------
data = CSV.read(DATA, DataFrame)
data.lambda = data.g

# --- simulation ----------------------------------------------------------------
p = Params()
runs = DataFrame[]
for (i, ry) in enumerate(0.0001:0.2:5)
    #     x    y    z                         r_x    r_y r_z                       d    λ
    u0 = [0.1, 0.2, rand(Normal(0.7, 0.009)), 0.001, ry, rand(Normal(0.2, 0.01)), 0.0, 1.0]
    df = simulate(p, u0; tmax=6.0)
    df.lineage .= i
    df.arr .= df.lambda[end] < 0.1 * p.lambda0
    push!(runs, df)
end
sim = vcat(runs...)
for g in (:x, :y, :z)
    sim[!, g] ./= mean(sim[sim.arr, g])
end

# --- plot ----------------------------------------------------------------------
# arrested lineages in the (y, v) plane: + start, dot end (`last` picks the end)
function plot_arrested(ax, df, v, last)
    for d in groupby(df[df.arr .== 1, :], :lineage)
        ax.plot(d.y, d[!, v], "-", color=ARRESTED, alpha=0.4, lw=1)
        ax.plot([d.y[1]], [d[1, v]], "k+")
        ax.plot(last(d.y), last(d[!, v]), "k.")
    end
end
endpoint(v) = [v[end]]
timeavg(v)  = [mean(v[2:end])]

fig, axs = subplots(figsize=(5.5, 6.6), nrows=3, ncols=2)

s = summarise(data)
axs[0, 0].plot(s.yv, s.lambda, "k.")
plot_arrested(axs[1, 0], data[data.time .< 50, :], :x, endpoint)
plot_arrested(axs[2, 0], data[data.time .< 50, :], :z, timeavg)
axs[0, 0].set_xlabel("Initial TetR velocity (a.u.)")
axs[1, 0].set_xlabel("TetR (a.u.)")
axs[2, 0].set_xlabel("TetR (a.u.)")
axs[0, 0].set_ylabel("Final growth rate")
axs[1, 0].set_ylabel("TetA (a.u.)")
axs[2, 0].set_ylabel("Constitutive (a.u.)")
axs[0, 0].set_title("Experiment")

s = summarise(sim)
axs[0, 1].plot(s.yv, s.lambda, "k.")
plot_arrested(axs[1, 1], sim, :x, endpoint)
plot_arrested(axs[2, 1], sim, :z, endpoint)
axs[0, 1].set_ylim(-0.05, 1.05)
axs[0, 1].set_xlabel("Initial velocity \$\\dot{y}(0)\$")
axs[1, 1].set_xlabel("\$y\$")
axs[2, 1].set_xlabel("\$y\$")
axs[0, 1].set_ylabel("Final \$\\lambda\$")
axs[1, 1].set_ylabel("\$x\$")
axs[2, 1].set_ylabel("\$z\$")
axs[0, 1].set_title("Model")

for j in 1:2, i in 1:3
    panel_label(axs[i-1, j-1], "ABCDEF"[3(j-1)+i:3(j-1)+i])
end

save(fig, "fig5.pdf")
