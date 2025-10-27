# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

"""
#:: FUNCTIONS TO RUN PHOTOLYSIS LOGIC DEPEND ON REACTION TYPE
#-------------------------------------------------------------------------------------------

# SUNFUNCTIONS:
- SPD_logic_multiple: Simple photodissociation logic
- DPD_logic_multiple: Double photodissociation logic
- SPI_logic_multiple: Simple photoionisation logic
- DPI_logic_multiple: Double photoionisation logic
- DiPI_logic_multiple: Dissociative photoionisation logic
"""

function SPD_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::Real, sun_tuple)
    product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), species_name)
    reaction = PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, true)
    final_speeds_light, final_speeds_heavy = multiple_photodissociation(reaction, energy_vector)
    
    product_types = [product_names[2], product_names[3]]
    product_velocities = [final_speeds_heavy, final_speeds_light]

    return product_velocities, product_types
end

function DPD_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], false)
    final_speeds_light_1, final_speeds_heavy_old = multiple_photodissociation(reaction, energy_vector)

    reaction = PhotoReaction(energy_to_J(tsh_energy[2]), mean(Float32[norm(p) for p in final_speeds_heavy_old]), sun_tuple, species_name[2], false)
    final_speeds_light_2, final_speeds_heavy = multiple_photodissociation(reaction, energy_vector)

    final_speeds_light_2_norms = [norm(p) for p in final_speeds_light_2]
    final_speeds_light_1_norms = [norm(p) for p in final_speeds_light_1]
    final_speeds_heavy_norms = [norm(p) for p in final_speeds_heavy]

    data_speeds = DataFrame(
    Product = ["O", "AVG H"],
    Mean_Speed = [mean(final_speeds_heavy_norms), mean([mean(final_speeds_light_1_norms), mean(final_speeds_light_2_norms)])],
    Median_Speed = [median(final_speeds_heavy_norms), mean([median(final_speeds_light_1_norms), median(final_speeds_light_2_norms)])],
    STD_half = [std(final_speeds_heavy_norms)/2, maximum([std(final_speeds_light_1_norms)/2, std(final_speeds_light_2_norms)/2])]
    )

    println(data_speeds)
    
    product_types = ["O", "H", "H"]
    product_velocities = [final_speeds_heavy, final_speeds_light_1, final_speeds_light_2]

    return product_velocities, product_types
end

function SPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, true)
    product_velocities = [multiple_photoionisation(reaction, energy_vector)]
    product_types = [species_name[2]]

    return product_velocities, product_types
end

function DPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], true)
    final_speeds_ion = multiple_photoionisation(reaction, energy_vector)

    reaction = PhotoReaction(energy_to_J(tsh_energy[2]), Float32(mean([norm(p) for p in final_speeds_ion])), sun_tuple, species_name[2], true)
    product_velocities = [multiple_photoionisation(reaction, energy_vector)]
    product_types = [species_name[2][2]]

    return product_velocities, product_types
end

function DiPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], true)
    final_speeds_ion = multiple_photoionisation(reaction, energy_vector)

    product_names = map(s -> replace(s, r"\(\+\)" => ""), species_name[2])
    
    reaction = PhotoReaction(energy_to_J(tsh_energy[2]), Float32(mean([norm(p) for p in final_speeds_ion])), sun_tuple, product_names, true)
    final_speeds_light, final_speeds_heavy = multiple_photodissociation(reaction, energy_vector)
    
    product_velocities = [final_speeds_heavy, final_speeds_light]
    product_types = [species_name[2][2], species_name[2][3]]

    return product_velocities, product_types
end


"""
#:: FUNCTION: call_photodestruction_logic_multiple(current_reaction, photochemical_info, parent_velocity, sun_tuple, photon_energy)
#-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Call the apropriate photodestruction logic providing the relevant photochemical information

# INPUTS:
- current_reaction
- photochemical_info
- parent_velocity: 3D tuple with aprent velocity
- sun_tuple: 3D tuple containing direction of incident photons
- photon_energy: in J
"""
 
function call_photodestruction_logic_multiple(current_reaction::CurrentReaction, photochemical_info::Species, parent_velocity::Float32, sun_tuple::Tuple{Float32, Float32, Float32}, energy_vector::Vector{Float32})

    pr = current_reaction.present_reaction
    r_idx = current_reaction.reaction_index

    tsh_energies = get_tsh_energies(photochemical_info)
    species_names = get_species_names(photochemical_info)


    # Simple photodissociations
    if pr == "SPD"
        return SPD_logic_multiple(energy_vector, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)
    
    # H2O Double photodissociation
    elseif pr == "DPD" 
        return DPD_logic_multiple(energy_vector, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)
    
    # Simple photoionisations
    elseif pr == "SPI"
        return SPI_logic_multiple(energy_vector, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)

    # Double electron ejection of negative atomic hydrogen, equivalente to a double photoionization process
    elseif pr == "DPI"
        return DPI_logic_multiple(energy_vector, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)

    # Dissociative photoionisations
    elseif pr == "DiPI"
        return DiPI_logic_multiple(energy_vector, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)
    end

end

function call_photodestruction_logic_multiple(current_reaction::CurrentReaction, photochemical_info::Species, parent_velocity::Real, sun_tuple::Nothing, energy_vector::AbstractVector{<:Real})
    return call_photodestruction_logic_multiple(current_reaction, photochemical_info, Float32(parent_velocity), random_unit_tuple(), Float32.(energy_vector))
end

function call_photodestruction_logic_multiple(current_reaction::CurrentReaction, photochemical_info::Species, parent_velocity::Real, sun_tuple::AbstractVector{<:Real}, energy_vector::AbstractVector{<:Real})
    st = Tuple(Float32.(sun_tuple))
    return call_photodestruction_logic_multiple(current_reaction, photochemical_info, Float32(parent_velocity), st, Float32.(energy_vector))
end


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
export photodestruction_logic_multiple