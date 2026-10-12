# Fig. 1A: growth rate vs TetR for each experimental lineage.
include(joinpath(@__DIR__, "..", "src", "common.jl"))

data = CSV.read(DATA, DataFrame)

fig, ax = subplots(figsize=(1.5, 1.5))
for (arr, col) in ((1, "C1"), (0, "C0"))
    for lin in unique(data[data.arr .== arr, :lineage])
        d = data[(data.lineage .== lin) .& (data.time .< 50), :]
        ax.plot(d.y, d.g, "-", color=col, alpha=0.3, lw=0.4)
        ax.plot([d.y[1]], [d.g[1]], "+", color=col)
        ax.plot([d.y[end]], [d.g[end]], ".", color=col)
    end
end
ax.plot([0, 2000], [0, 0], "k--")
ax.set_yticks([0, 1])
ax.set_ylabel(L"$\lambda/\langle\lambda_0\rangle$")
ax.set_xlabel(L"TetR \,(a.u.)")
save(fig, "fig1a.pdf")
