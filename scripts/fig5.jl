# Fig. 5: experimental lineages (left) against simulations over a range of
# initial TetR resources r_y(0) (right).
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
# arrested lineages in the (y, v) plane; `last` picks the end marker
function plot_arrested(ax, df, v, last)
    for d in groupby(df[df.arr .== 1, :], :lineage)
        ax.plot(d.y, d[!, v], "-", color="C1", alpha=0.4)
        ax.plot([d.y[1]], [d[1, v]], "k+")
        ax.plot(last(d.y), last(d[!, v]), "k.")
    end
end
endpoint(v) = [v[end]]
timeavg(v)  = [mean(v[2:end])]

fig, axs = subplots(figsize=(6, 7), nrows=3, ncols=2)

s = summarise(data)
axs[0, 0].plot(s.yv, s.lambda, "k.")
axs[0, 0].set_xlabel(L"y'(0)")
axs[0, 0].set_ylabel(L"\lambda(t_f)")

early = data[data.time .< 50, :]
plot_arrested(axs[1, 0], early, :x, endpoint)
axs[1, 0].set_xlabel(L"TetR")
axs[1, 0].set_ylabel(L"TetA")
plot_arrested(axs[2, 0], early, :z, timeavg)
axs[2, 0].set_xlabel(L"TetR")
axs[2, 0].set_ylabel("constitutively expr. gene")

s = summarise(sim)
axs[0, 1].plot(s.yv, s.lambda, "k.")
axs[0, 1].set_ylim(0, 1)
axs[0, 1].set_xlabel(L"r_y^0")

plot_arrested(axs[1, 1], sim, :x, endpoint)
axs[1, 1].set_xlabel(L"y")
axs[1, 1].set_ylabel(L"x")
plot_arrested(axs[2, 1], sim, :z, endpoint)
axs[2, 1].set_xlabel(L"y")
axs[2, 1].set_ylabel(L"z")

for i in 1:3, (j, col) in enumerate("AB")
    axs[i-1, j-1].text(0.7, 0.9, "($col.$i)", transform=axs[i-1, j-1].transAxes,
                       fontsize=15, va="top", ha="left")
end
fig.tight_layout()
save(fig, "fig5.pdf")
