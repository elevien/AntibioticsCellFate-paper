# Resource-allocation model. State is u = [x..., r..., d, λ]: protein levels,
# their resources, internal drug and growth rate. The first gene is the efflux pump.

@kwdef struct Params
    sectors::Vector{Symbol} = [:P, :P, :Q]   # x (TetA), y (TetR), z (constitutive)
    A::Matrix{Float64} = [1.0 -1.0 0.0;
                          0.0  1.0 0.0;
                          0.0  0.0 1.0]
    lambda0::Float64 = 1.0
    Kd::Float64      = 666.67
    h::Float64       = 5.0
    J::Float64       = 1.0
    beta::Float64    = 5333.33
    Kx::Float64      = 666.67
    D::Float64       = 2666.67
    gamma::Float64   = 1.0
    # λ relaxes to growth(d) at this rate; it stands in for the algebraic
    # relation λ = λ(d). Euler needs dt * relax < 2.
    relax::Float64   = 100.0
end

growth(d, p) = p.lambda0 / (1 + (d / p.Kd)^p.h)

function rhs!(du, u, p)
    n = length(p.sectors)
    x, r = u[1:n], u[n+1:2n]
    d, λ = u[2n+1], u[2n+2]
    # P-sector resources track the growth rate, Q-sector resources do not
    target = [s == :P ? λ / p.lambda0 : 1.0 for s in p.sectors]

    du[1:n]    .= λ .* (r .- x)
    du[n+1:2n] .= λ^p.gamma .* (p.A * (target .- r))
    du[2n+1]    = p.J * (p.D - d) - (p.beta * x[1] / (p.Kx + x[1]) + λ) * d
    du[2n+2]    = p.relax * (growth(d, p) - λ)
end

"""
    simulate(p, u0; tmax, dt=0.01)

Fixed-step Euler integration. Returns a DataFrame with columns
`time, x, y, z, r_x, r_y, r_z, d, lambda` (fewer genes give fewer columns).
"""
function simulate(p::Params, u0; tmax, dt=0.01)
    n = length(p.sectors)
    @assert length(u0) == 2n + 2
    times = collect(0:dt:tmax)
    out = Matrix{Float64}(undef, length(times), 2n + 2)

    u, du = float.(u0), zeros(2n + 2)
    out[1, :] .= u
    for k in 2:length(times)
        rhs!(du, u, p)
        u .+= dt .* du
        out[k, :] .= u
    end

    genes = ["x", "y", "z"][1:n]
    df = DataFrame(out, [genes; "r_" .* genes; "d"; "lambda"])
    insertcols!(df, 1, :time => times)
    return df
end
