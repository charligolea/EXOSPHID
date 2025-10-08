using BenchmarkTools
using Distributions

include("../src/ExosPhID.jl")
using .ExosPhID

reaction_occurence, reaction_name, product_types, product_velocities = true, "", [], []

dt = 1.0f10 # s, Large dt to make sure interaction happens
solar_activity = 0.0f0
velocities = Dict("H2O" => 590.0f0, "OH" => 605.0f0, "H2" => 1750.0f0, "H" => 2500.0f0, "H(-)"=> 2500.0f0, "HO2" => 425.0f0, "H2O2" => 435.0f0, "He" => 1250.0f0, "Ne" => 560.0f0)

for parent_name in ["H2O"] # Modify at will

    print("\nStarting Numerical Benchmark for $(parent_name)\n")

    parent_vector = fill(parent_name, 1)

    for (index, parent_type) in enumerate(parent_vector)

        if parent_type in ("H2O", "H2", "OH", "H", "H(-)", "HO2", "H2O2", "He", "Ne")
            θ_parent, φ_parent = 2 * π * rand(), π * rand()
            parent_velocity =  velocities[parent_type] .*  (Float32(sin(φ_parent)) * Float32(cos(θ_parent)), Float32(sin(φ_parent)) * Float32(sin(θ_parent)), Float32(cos(φ_parent)))
        else
            print("Parent species is not valid")
        end

        bm = @benchmark photodestruction($solar_activity, $dt, $parent_type, $parent_velocity, nothing)
        display(bm)
    end

    print("\nNumerical Benchmark completed for $(parent_name)\n")

end