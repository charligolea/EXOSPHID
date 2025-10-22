"""include("../database/photodatabase.jl")
using .photodatabase"""

include("../database/species/general_construct.jl")

"""
#:: Various functions to unpack photochemical_info
#-------------------------------------------------------------------------------------------

# SUNFUNCTIONS:
- get_tsh_energies(photo_info) = unpack dissociation / ionisation energies for species of interest
- get_species_names(photo_info) = unpack parent and product species names for reactions for species of interest
- unpack_present_reaction(reaction_info) = unpack dissociation / ionisation energies for reaction of interest
- unpack_reaction_index(reaction_info) = unpack parent and product species names for reaction of interest
- energy_to_J(energy) = converte enrgy to J
"""

get_tsh_energies(photo_info) = photo_info.tsh_energies  # eV
get_species_names(photo_info) = photo_info.species_names
energy_to_J(energy) = 1.602f-19 * energy


"""
#:: FUNCTIONS TO RUN PHOTOLYSIS LOGIC DEPEND ON REACTION TYPE
#-------------------------------------------------------------------------------------------

# SUNFUNCTIONS:
- SPD_logic: Simple photodissociation logic
- DPD_logic: Double photodissociation logic
- SPI_logic: Simple photoionisation logic
- DPI_logic: Double photoionisation logic
- DiPI_logic: Dissociative photoionisation logic
"""

function SPD_logic(photon_energy::Float32, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), species_name)
    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, false)
    final_speeds_light, final_speeds_heavy = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)
    
    product_types = [product_names[2], product_names[3]]
    product_velocities = [final_speeds_heavy, final_speeds_light]

    return product_velocities, product_types
end

function DPD_logic(photon_energy::Float32, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], false)
    final_speeds_light_1, final_speeds_heavy_old = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)

    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy[2]), map(Float32, final_speeds_heavy_old), sun_tuple, species_name[2], false)
    final_speeds_light_2, final_speeds_heavy = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)
    
    product_types = ["O", "H", "H"]
    product_velocities = [final_speeds_heavy, final_speeds_light_1, final_speeds_light_2]

    return product_velocities, product_types
end

function SPI_logic(photon_energy::Float32, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, false)
    product_velocities = [SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)]
    product_types = [species_name[2]]

    return product_velocities, product_types
end

function DPI_logic(photon_energy::Float32, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], false)
    final_speeds_ion = SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)

    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy[2]), map(Float32, final_speeds_ion), sun_tuple, species_name[2], false)
    product_velocities = [SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)]
    product_types = [species_name[2][2]]

    return product_velocities, product_types
end

function DiPI_logic(photon_energy::Float32, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], false)
    final_speeds_ion = SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)

    product_names = map(s -> replace(s, r"\(\+\)" => ""), species_name[2])
    
    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy[2]), map(Float32, final_speeds_ion), sun_tuple, product_names, false)
    final_speeds_light, final_speeds_heavy = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)
    
    product_velocities = [final_speeds_heavy, final_speeds_light]
    product_types = [species_name[2][2], species_name[2][3]]

    return product_velocities, product_types
end
 

"""
#:: FUNCTION: call_photodestruction_logic(current_reaction, photochemical_info, parent_velocity, sun_tuple, photon_energy)
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

function call_photodestruction_logic(current_reaction, photochemical_info, parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32}, photon_energy::Real)

    pr = current_reaction.present_reaction
    r_idx = current_reaction.reaction_index

    tsh_energies = get_tsh_energies(photochemical_info)
    species_names = get_species_names(photochemical_info)

    # Simple photodissociations
    if pr == "SPD"
        return SPD_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)
    
    # Double photodissociation
    elseif pr == "DPD" 
        return DPD_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)
    
    # Simple photoionisations
    elseif pr == "SPI"
        return SPI_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)

    # Double electron ejection of negative atomic hydrogen, equivalent to a double photoionization process
    elseif pr == "DPI"
        return DPI_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)

    # Dissociative photoionisations
    elseif pr == "DiPI"
        return DiPI_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], parent_velocity, sun_tuple)
    end

end