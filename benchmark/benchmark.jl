module benchmark

using BenchmarkTools
using Distributions
using Profile

include("../src/ExosPhID.jl")
using .ExosPhID

include("../src/database/photodatabase.jl")
using .photodatabase

export photobenchmark
export allocations

"""
#:: FUNCTION: photobenchmark(dt, solar_activity, parent_name)

# OBJECTIVE: Evaluate computational performance of ExosPhID

# INPUTS:
- dt: Time Window in seconds
- solar_activity: From 0 (Quiet Sun) to 1 (Active Sun)
- parent_name: String with parent molecule type. See exosphid_species variable in photodatabase.jl for possible species
"""

function photobenchmark(dt::Float32, solar_activity::Float32, parent_name::String)

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


"""
#:: FUNCTION: allocations(dt, solar_activity, parent_name)

# OBJECTIVE: Evaluate allocations in detail

# INPUTS:
- dt: Time Window in seconds
- solar_activity: From 0 (Quiet Sun) to 1 (Active Sun)
- parent_name: String with parent molecule type. See exosphid_species variable in photodatabase.jl for possible species
"""

function allocations(dt::Float32, solar_activity::Float32, parent_name::String)

    println("\nStarting Allocation Profiler for $(parent_name)\n")

    velocities = Dict("H2O" => 590, "OH" => 605, "H2" => 1750, "H" => 2500, 
                      "H(-)" => 2500, "HO2" => 425, "H2O2" => 435, "He" => 1250, "Ne" => 560)
    parent_velocity = velocities[parent_name]

    Profile.clear()
    @profile photodestruction(solar_activity, dt, parent_name, parent_velocity, nothing)

    println("\nAllocation hotspots:")
    Profile.print(format=:flat, sortedby=:allocs, maxdepth=20)
end

function allocations(dt::Real, solar_activity::Real, parent_name::String)
    return allocations(Float32(dt), Float32(solar_activity), parent_name)
end

end