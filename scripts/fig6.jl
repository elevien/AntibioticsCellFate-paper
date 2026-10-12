# Fig. 6: critical initial resource as a function of external drug D.
include(joinpath(@__DIR__, "..", "src", "common.jl"))

#     x       y    z    r_x    r_y  r_z  d    λ
u0 = [0.0001, 0.3, 1.0, 0.001, 0.6, 0.6, 0.0, 1.0]
r_range = 0.0001:0.1:10

# NOTE: the swept initial condition is r_x(0); the axis label follows the paper.
arrested(D, r) =
    simulate(Params(D=D), [u0[1:3]; r; u0[5:8]]; tmax=5.0).lambda[end] < 0.05

# largest initial resource that still ends up arrested (nothing if all recover)
function critical_r(D)
    i = findlast(r -> arrested(D, r), r_range)
    return isnothing(i) ? nothing : r_range[i]
end

Ds, r_crit = Float64[], Float64[]
for D in range(1600.0, 14667.0, length=20)
    rc = critical_r(D)
    isnothing(rc) && continue
    push!(Ds, D)
    push!(r_crit, rc)
end

fig, ax = subplots(figsize=(3, 3))
ax.plot(Ds, r_crit, ".")
ax.loglog(Ds, 1.5 * r_crit[end] .* (Ds ./ Ds[end]) .^ 2, "k--", label=L"$D^2$")
ax.legend(frameon=false)
ax.set_xticks([1e3, 1e4])
ax.set_ylim(1e-1, 20)
ax.set_xlabel(L"D")
ax.set_ylabel(L"r_y(0)^*")
save(fig, "fig6.pdf")
