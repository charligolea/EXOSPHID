module validation

export validate_exothermic_velocities

include("../src/photoreactions/SimplePhotodissociation.jl")
using .SimplePhotodissociation
include("../src/photoreactions/SimplePhotoionisation.jl")
using .SimplePhotoionisation

include("../src/database/solar_spectrum.jl")
using .solar_spectrum

include("../src/database/photodatabase.jl")
include("../src/photoreactions/photodestruction_logic.jl")
include("../src/photoreactions/photodestruction_logic_multiple.jl")


""" 
This script simulates every possible reaction for a given species a given number of times (we assume that a reaction is always happening, so infinite dt!).
The idea is to run it many times (thousands to millions) to compare to reported values of exothermic velocities in literature.

# Arguments
- sun_mode = From 0 (Quiet Sun) to 1 (Active Sun)
- parent_type:: String with parent molecule type ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2")
- num_reactions:: Number of photon interactions to simulate for every possible reaction    
"""

function validate_exothermic_velocities(parent_type::String, sun_mode::Float32, num_reactions::Int)

    @assert 0.0 <= sun_mode <= 1.0 "Solar activity must be in (0,1)!"
    @assert parent_type in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne") "Invalid parent species: $parent_type"

    velocities = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, "He" => 1250.0, "Ne" => 560.0)
    parent_velocity =  velocities[parent_type]

    photochemical_info = get_species_photochemical_info(parent_type)
    wavelength_range = get_wavelength_range(photochemical_info)
    reaction_names = photochemical_info[4]
    reaction_probabilities = photochemical_info[3]
    reaction_types = photochemical_info[5]

    for (w_index, wvl_range) in enumerate(wavelength_range)

        print("Wavelength Range: ", wvl_range, "\n\n")

        wvl_vector, energy_vector = flux_outputs(parent_type, wvl_range, sun_mode, num_reactions)

        for (r_index, rn) in enumerate(reaction_names)
            if reaction_probabilities[r_index][w_index] > 0.0
                print("Simulating reaction: ", reaction_names[r_index], "\n")
                photoreaction_characteristics = (reaction_types[r_index], "", r_index, "")
                product_velocities, product_types = call_photodestruction_logic_multiple(photoreaction_characteristics, photochemical_info, parent_velocity, nothing, energy_vector)
            end
        end
        print("\n\n")

    end

end

function validate_exothermic_velocities(parent_type::String, sun_mode::Real, num_reactions::Int)
    return validate_exothermic_velocities(parent_type, Float32(sun_mode), num_reactions)
end

end