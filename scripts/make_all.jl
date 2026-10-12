# Run every figure script, each in its own module.
for n in 1:6
    @eval module $(Symbol(:Fig, n))
        include(joinpath(@__DIR__, $("fig$n.jl")))
    end
end
