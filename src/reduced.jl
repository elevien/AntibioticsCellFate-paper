# Nullclines, fixed points and linear stability of the reduced (x, d) and
# (x, r, d) models. Fixed points are parametrised by the internal drug level d.

growth_prime(d, p) =
    d <= 0 ? 0.0 : -p.lambda0 * p.h / p.Kd * (d / p.Kd)^(p.h - 1) / (1 + (d / p.Kd)^p.h)^2

# x-nullcline
x_nc(d, p) = growth(max(0.0, d), p) / p.lambda0

# d-nullcline, solved for x
function d_nc(d, D, p)
    a = p.J * (D / d - 1) - growth(d, p)
    return a * p.Kx / (p.beta - a)
end

# external drug level D whose fixed point sits at internal level d
function D_at_fixed_point(d, p)
    d <= 0 && return 0.0
    x = x_nc(d, p)
    return d * (1 + (p.beta * x / (p.Kx + x) + growth(d, p)) / p.J)
end

# Jacobian entries at the fixed point
function jac_terms(d, p)
    x = x_nc(d, p)
    λ = growth(d, p)
    Λ = -growth_prime(d, p)
    C = p.J + p.beta * x / (p.Kx + x) + λ - Λ * d
    μ = p.beta * p.Kx * d / (p.Kx + x)^2
    return λ, Λ, C, μ
end

function max_re_2d(d, p)
    d <= 0 && return -growth(0.0, p)
    λ, Λ, C, μ = jac_terms(d, p)
    return maximum(real, eigvals([-λ -λ*Λ; -μ -C]))
end

function max_re_3d(d, γ, p)
    d <= 0 && return -growth(0.0, p)
    λ, Λ, C, μ = jac_terms(d, p)
    B = λ^γ
    return maximum(real, eigvals([-λ  λ    0.0;
                                  0.0 -B  -B*Λ;
                                  -μ  0.0  -C]))
end

stability_2d(d, p) = max_re_2d(d, p) < 0 ? :stable : :saddle

# fixed points (d, x, stability) of the 2D model at external drug D,
# located by sign changes of D_at_fixed_point over the grid d_grid
function fixed_points_2d(D, d_grid, p)
    fps = Tuple{Float64,Float64,Symbol}[]
    Δ = [D_at_fixed_point(d, p) for d in d_grid] .- D
    for i in 1:length(d_grid)-1
        if Δ[i] * Δ[i+1] < 0
            d = d_grid[i] + Δ[i] / (Δ[i] - Δ[i+1]) * (d_grid[i+1] - d_grid[i])
            push!(fps, (d, x_nc(d, p), stability_2d(d, p)))
        end
    end
    return fps
end

# Euler integration of the 2D model; returns (t, λ, x)
function simulate_2d(D, p; x0=0.01, d0=0.0, tmax=30.0, dt=0.01)
    ts = collect(0:dt:tmax)
    λs, xs = similar(ts), similar(ts)
    x, d = x0, d0
    for k in eachindex(ts)
        λ = growth(d, p)
        λs[k], xs[k] = λ, x
        x, d = x + dt * λ * (x_nc(d, p) - x),
               d + dt * (p.J * (D - d) - (p.beta * x / (p.Kx + x) + λ) * d)
    end
    return ts, λs, xs
end
