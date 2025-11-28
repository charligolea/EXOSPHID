using BenchmarkTools
using Profile

if !isdefined(Main, :EXOSPHID)
    include(joinpath(@__DIR__, "..", "src", "EXOSPHID.jl"))
    using .EXOSPHID
end


"""
    photobenchmark(dt, solar_activity, parent_name)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Evaluate computational performance of EXOSPHID

# INPUTS:
- `dt::Float32` -> Time Window in seconds
- `solar_activity::Float32` -> From 0 (Quiet Sun) to 1 (Active Sun)
- `parent_name::String` -> String with parent molecule type. See exosphid_species variable 
    in photodatabase.jl for possible species
"""
function photobenchmark(dt::Float32, solar_activity::Float32, parent_name::String)

    print("\nStarting Numerical Benchmark for $(parent_name)\n")
    parent_velocity =  EXOSPHID.velocities[parent_name]

    bm = @benchmark photodestruction($solar_activity, $dt, $parent_name, 
                                    $parent_velocity, nothing)
    display(bm)

    print("\nNumerical Benchmark completed for $(parent_name)\n")

end

function photobenchmark(dt::Real, solar_activity::Real, parent_name::String)
    return photobenchmark(Float32(dt), Float32(solar_activity), parent_name)
end


"""
    allocations(dt, solar_activity, parent_name)
-------------------------------------------------------------------------------------------

# OBJECTIVE: Evaluate allocations in detail

# INPUTS:
- `dt::Float32` -> Time Window in seconds
- `solar_activity::Float32` -> From 0 (Quiet Sun) to 1 (Active Sun)
- `parent_name::String` -> String with parent molecule type. See exosphid_species variable 
    in photodatabase.jl for possible species
"""
function allocations(dt::Float32, solar_activity::Float32, parent_name::String)

    println("\nStarting Allocation Profiler for $(parent_name)\n")
    parent_velocity = EXOSPHID.velocities[parent_name]

    Profile.clear()
    @profile photodestruction(solar_activity, dt, parent_name, parent_velocity, nothing)

    println("\nAllocation hotspots:")
    Profile.print(format=:flat, sortedby=:allocs, maxdepth=20)
end

function allocations(dt::Real, solar_activity::Real, parent_name::String)
    return allocations(Float32(dt), Float32(solar_activity), parent_name)
end

export photobenchmark
export allocations