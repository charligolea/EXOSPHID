"""
    const velocities

- Sample values of parent velocities for the different EXOSPHID species.
- Determined from thermal velocity distributions for an assumed temperature of 250 K 
    (average for the Moon)
"""
const velocities = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, 
                        "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, 
                        "He" => 1250.0, "Ne" => 560.0, "Ar" => 395.0)

""" 
    validate_exothermic_velocities(parent_type, solar_activity, num_reactions)
-------------------------------------------------------------------------------------------

This script simulates every possible reaction for a given species a given number of times 
    (we assume that a reaction is always happening, so infinite dt!).
The idea is to run it many times (thousands to millions) to compare to reported values of 
    exothermic velocities in literature.

# Arguments
- `solar_activity::Float32` -> From 0 (Quiet Sun) to 1 (Active Sun)
- `parent_type::String` -> String with parent molecule type in `exosphid_species`
- `num_reactions::Int` -> Number of photon interactions to simulate for every possible 
    reaction    
"""
function validate_exothermic_velocities(parent_type::String, solar_activity::Float32, 
                                        num_reactions::Int)

    parent_velocity =  velocities[parent_type]

    photochemical_info = get_species_photochemical_info(parent_type)
    wrs = get_wavelength_range(photochemical_info)
    rns = photochemical_info.reaction_names
    rps = photochemical_info.reaction_probabilities
    rts = photochemical_info.reaction_types

    for (w_index, wvl_range) in enumerate(wrs)

        print("Wavelength Range: ", wvl_range, "\n\n")

        wvl_vector, energy_vector = 
            flux_outputs(parent_type, wvl_range, solar_activity, num_reactions)

        for (r_index, rn) in enumerate(rns)
            if rps[r_index][w_index] > 0.0
                print("Simulating reaction: ", rns[r_index], "\n")
                current_reaction = CurrentReaction(rts[r_index], "", r_index, wvl_range)
                product_velocities, product_types = 
                    call_photodestruction_logic_multiple(current_reaction, 
                    photochemical_info, parent_velocity, nothing, energy_vector)
            end
        end
        print("\n\n")

    end

end

function validate_exothermic_velocities(parent_type::String, solar_activity::Real, 
                                        num_reactions::Int)
    return validate_exothermic_velocities(parent_type, Float32(solar_activity), 
        num_reactions)
end

export validate_exothermic_velocities