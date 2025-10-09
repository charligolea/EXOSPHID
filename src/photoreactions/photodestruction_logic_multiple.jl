using DataFrames
using Statistics 
using LinearAlgebra

get_tsh_energies(photo_info) = photo_info[1]  # Already converted to J
get_species_names(photo_info) = photo_info[2]
unpack_present_reaction(reaction_info) = reaction_info[1]
unpack_reaction_index(reaction_info) = reaction_info[3]
energy_to_J(energy) = 1.602f-19 * energy

function SPD_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::Real, sun_tuple)
    product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), species_name)
    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, true)
    final_speeds_light, final_speeds_heavy = SimplePhotodissociation.multiple_photodissociation(reaction, energy_vector)
    
    product_types = [product_names[2], product_names[3]]
    product_velocities = [final_speeds_heavy, final_speeds_light]

    return product_velocities, product_types
end

function DPD_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], false)
    final_speeds_light_1, final_speeds_heavy_old = SimplePhotodissociation.multiple_photodissociation(reaction, energy_vector)

    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy[2]), mean(Float32[norm(p) for p in final_speeds_heavy_old]), sun_tuple, species_name[2], false)
    final_speeds_light_2, final_speeds_heavy = SimplePhotodissociation.multiple_photodissociation(reaction, energy_vector)

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
    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, true)
    product_velocities = [SimplePhotoionisation.multiple_photoionisation(reaction, energy_vector)]
    product_types = [species_name[2]]

    return product_velocities, product_types
end

function DPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], true)
    final_speeds_ion = SimplePhotoionisation.multiple_photoionisation(reaction, energy_vector)

    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy[2]), Float32(mean([norm(p) for p in final_speeds_ion])), sun_tuple, species_name[2], true)
    product_velocities = [SimplePhotoionisation.multiple_photoionisation(reaction, energy_vector)]
    product_types = [species_name[2][2]]

    return product_velocities, product_types
end

function DiPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = SimplePhotoionisation.PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], true)
    final_speeds_ion = SimplePhotoionisation.multiple_photoionisation(reaction, energy_vector)

    product_names = map(s -> replace(s, r"\(\+\)" => ""), species_name[2])
    
    reaction = SimplePhotodissociation.PhotoReaction(energy_to_J(tsh_energy[2]), Float32(mean([norm(p) for p in final_speeds_ion])), sun_tuple, product_names, true)
    final_speeds_light, final_speeds_heavy = SimplePhotodissociation.multiple_photodissociation(reaction, energy_vector)
    
    product_velocities = [final_speeds_heavy, final_speeds_light]
    product_types = [species_name[2][2], species_name[2][3]]

    return product_velocities, product_types
end
 
function call_photodestruction_logic_multiple(photoreaction_characteristics, photochemical_info, parent_velocity, sun_tuple, energy_vector)

    present_reaction = unpack_present_reaction(photoreaction_characteristics)
    reaction_index = unpack_reaction_index(photoreaction_characteristics)

    tsh_energies = get_tsh_energies(photochemical_info)
    species_names = get_species_names(photochemical_info)


    # Simple photodissociations
    if present_reaction == "SPD"
        return SPD_logic_multiple(energy_vector, species_names[reaction_index], tsh_energies[reaction_index], parent_velocity, sun_tuple)
    
        # H2O Double photodissociation
    elseif present_reaction == "DPD" 
        return DPD_logic_multiple(energy_vector, species_names[reaction_index], tsh_energies[reaction_index], parent_velocity, sun_tuple)
    
    # Simple photoionisations
    elseif present_reaction == "SPI"
        return SPI_logic_multiple(energy_vector, species_names[reaction_index], tsh_energies[reaction_index], parent_velocity, sun_tuple)

    # Double electron ejection of negative atomic hydrogen, equivalente to a double photoionization process
    elseif present_reaction == "DPI"
        return DPI_logic_multiple(energy_vector, species_names[reaction_index], tsh_energies[reaction_index], parent_velocity, sun_tuple)

    # Dissociative photoionisations
    elseif present_reaction == "DiPI"
        return DiPI_logic_multiple(energy_vector, species_names[reaction_index], tsh_energies[reaction_index], parent_velocity, sun_tuple)
    end

end