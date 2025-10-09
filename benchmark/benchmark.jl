module benchmark

using BenchmarkTools
using Distributions

include("../src/ExosPhID.jl")
using .ExosPhID

export photobenchmark

function photobenchmark(dt::Float32, solar_activity::Float32, parent_name::String)

    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"
    @assert 0.0 <= dt "dt must be positive!"
    @assert parent_name in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne") "Invalid parent species: $parent_name"

    reaction_occurence, reaction_name, product_types, product_velocities = true, "", [], []
    velocities = Dict("H2O" => 590, "OH" => 605, "H2" => 1750, "H" => 2500, "H(-)"=> 2500, "HO2" => 425, "H2O2" => 435, "He" => 1250, "Ne" => 560)

    print("\nStarting Numerical Benchmark for $(parent_name)\n")

    parent_velocity =  velocities[parent_name]

    bm = @benchmark photodestruction($solar_activity, $dt, $parent_name, $parent_velocity, nothing)
    display(bm)

    print("\nNumerical Benchmark completed for $(parent_name)\n")

end

function photobenchmark(dt::Real, solar_activity::Real, parent_name::String)
    return photobenchmark(Float32(dt), Float32(solar_activity), parent_name)
end

end