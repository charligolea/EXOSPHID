using JET
using Distributions

if !isdefined(Main, :EXOSPHID)
    include(joinpath(@__DIR__, "..", "src", "EXOSPHID.jl"))
    using .EXOSPHID
end

const velocities = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, "He" => 1250.0, "Ne" => 560.0)

dt = 1.0f10 # s, Large dt to make sure interaction happens
solar_activity = 0.0f0
velocities = Dict("H2O" => 590.0f0, "OH" => 605.0f0, "H2" => 1750.0f0, "H" => 2500.0f0, "H(-)"=> 2500.0f0, "HO2" => 425.0f0, "H2O2" => 435.0f0, "He" => 1250.0f0, "Ne" => 560.0f0)

parent_name = "H2O"
parent_velocity = velocities[parent_type]
parent_vector = fill(parent_name, 1)

print("\nStarting JET Test 1: Type Instability\n\n")

for (index, parent_type) in enumerate(parent_vector)
    @report_opt photodestruction(solar_activity, dt, parent_type, parent_velocity, nothing)    
end

print("End of JET Test 1: Type Instability\n\n\n")


print("Starting JET Test 2: Type Errors\n\n")

for (index, parent_type) in enumerate(parent_vector)
    @report_call photodestruction(solar_activity, dt, parent_type, parent_velocity, nothing)    
end

print("End of JET Test 2: Type Errors\n\n\n")

print("Starting JET Test 3: Full Analysis\n\n")

for (index, parent_type) in enumerate(parent_vector)
    @code_warntype photodestruction(solar_activity, dt, parent_type, parent_velocity, nothing)
end

print("End of JET Test 3: Full Analysis")