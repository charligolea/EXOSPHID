# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS TO RUN PHOTOLYSIS LOGIC DEPEND ON REACTION TYPE
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

"""
    SPD_logic_multiple(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple): 

# OBJECTIVE:
- Handles Simple PhotoDissociation (SPD) logic

# INPUTS:
- `energy_vector::Vector{Float32}` -> Energies of generated photons, in J
- `snames::NTuple{3, String}` -> 3D Tuple containing names of parent, and daughter (x2) species, POTENTIALLY INCLUDING ELECTRONIC STATES IN PARENTHESIS
- `tsh_energy::Float32` -> Bond/Dissociation Energy for parent species, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 2-element vector. Each element is a NTuple{3, Float64} -> Each 3D Tuple contains the velocity vector for the 2 daughter species specified in product_types 
- `product_types` -> 2-element vector containing String objects -> Name of the 2 daughter species. 
    Only considers "plain" atomic/molecular species, meaning no information about electronic states is given.
"""
function SPD_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::Real, sun_tuple)
    product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), species_name)
    reaction = PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, true)
    final_speeds_light, final_speeds_heavy = multiple_photodissociation(reaction, energy_vector)
    
    product_types = [product_names[2], product_names[3]]
    product_velocities = [final_speeds_heavy, final_speeds_light]

    return product_velocities, product_types
end


"""
    DPD_logic_multiple(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple): 

# OBJECTIVE:
- Handles Double PhotoDissociation (DPD) logic

# INPUTS:
- `energy_vector::Vector{Float32}` -> Energies of generated photons, in J
- `snames::NTuple{2, NTuple{3, String}}` -> Tuple congtaining 2 elements. 
    * Each of these elements is a 3D Tuple containing names of parent, and daughter (x2) species
    * The first Tuple will corrspond to the first dissociation, the other to the second consequent dissociation
    * POTENTIALLY INCLUDING ELECTRONIC STATES IN PARENTHESIS
- `tsh_energy::NTuple{2, Float32}` -> Bond/Dissociation Energy for first and for the total reaction, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 3-element vector. Each element is a NTuple{3, Float64} -> Each 3D Tuple contains the velocity vector for the 3 daughter species specified in product_types 
- `product_types` -> 3-element vector containing String objects -> Name of the 3 daughter species. 
    * Only considers "plain" atomic/molecular species, meaning no information about electronic states is given.
    * The first species will be the heavy product, the other 2 the lighter ones.
    * Currently, this reaction is only supported for the DPD of Water: H2O -> OH + H -> O + H + H
"""
function DPD_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], false)
    final_speeds_light_1, final_speeds_heavy_old = multiple_photodissociation(reaction, energy_vector)

    reaction = PhotoReaction(energy_to_J(tsh_energy[2]), mean(Float32[norm(p) for p in final_speeds_heavy_old]), sun_tuple, species_name[2], false)
    final_speeds_light_2, final_speeds_heavy = multiple_photodissociation(reaction, energy_vector)

    final_speeds_light_2_norms = [norm(p) for p in final_speeds_light_2]
    final_speeds_heavy_norms = [norm(p) for p in final_speeds_heavy]

    data_speeds = DataFrame(
    Product = ["O", "H"],
    Mean_Speed = [mean(final_speeds_heavy_norms), mean(final_speeds_light_2_norms)],
    Median_Speed = [median(final_speeds_heavy_norms), median(final_speeds_light_2_norms)],
    STD_half = [std(final_speeds_heavy_norms)/2, std(final_speeds_light_2_norms)/2]
    )

    println(data_speeds)
    
    product_types = ["O", "H", "H"]
    product_velocities = [final_speeds_heavy, final_speeds_light_1, final_speeds_light_2]

    return product_velocities, product_types
end


"""
    SPI_logic_multiple(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple): 

# OBJECTIVE:
- Handles Simple PhotoIonisation (SPI) logic

# INPUTS:
- `energy_vector::Vector{Float32}` -> Energies of generated photons, in J
- `snames::NTuple{3, String}` -> 3D Tuple containing names of parent (e.g H2O), ionised product (e.g. H2O(+)) and electron (e.g. e(-))
    * electronic states should not be included in this case
- `tsh_energy::Float32` -> Ionisation Energy in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 1-element vector, containing NTuple{3, Float64} -> The 3D Tuple contains the velocity vector for the ionised parent species
- `product_types` -> 1-element vector containing String object -> Name of the ionised parent
"""
function SPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{3, String}, tsh_energy::Float32, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy), parent_velocity, sun_tuple, species_name, true)
    product_velocities = [multiple_photoionisation(reaction, energy_vector)]
    product_types = [species_name[2]]

    return product_velocities, product_types
end


"""
    DPI_logic_multiple(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple): 

# OBJECTIVE:
- Handles Double PhotoIonisation (DPI) logic

# INPUTS:
- `energy_vector::Vector{Float32}` -> Energies of generated photons, in J
- `snames::NTuple{3, String}` -> Tuple containing 2 Tuple elements.
    * Both of these elements are 3D Tuples containing names of parent (e.g H2O), ionised product (e.g. H2O(+)) and electron (e.g. e(-))
    * Two tuples, for 2 ionisation reactions
    * electronic states should not be included in parenthesis in this case
- `tsh_energy::NTuple{2, Float32}` -> Ionisation Energy for first and for total reaction, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 1-element vector, containing NTuple{3, Float64} -> The 3D Tuple contains the velocity vector for the final ionised product
- `product_types` -> 1-element vector containing String object -> Name of the final ionised species
"""
function DPI_logic_multiple(energy_vector::Vector{Float32}, species_name::NTuple{2, NTuple{3, String}}, tsh_energy::NTuple{2, Float32}, parent_velocity::Real, sun_tuple)
    reaction = PhotoReaction(energy_to_J(tsh_energy[1]), parent_velocity, sun_tuple, species_name[1], true)
    final_speeds_ion = multiple_photoionisation(reaction, energy_vector)

    reaction = PhotoReaction(energy_to_J(tsh_energy[2]), Float32(mean([norm(p) for p in final_speeds_ion])), sun_tuple, species_name[2], true)
    product_velocities = [multiple_photoionisation(reaction, energy_vector)]
    product_types = [species_name[2][2]]

    return product_velocities, product_types
end


"""
    DiPI_logic_multiple(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple): 

# OBJECTIVE:
- Handles Dissociative PhotoIonisation (DiPI) logic
- The parent species is ionised, then dissociated.

# INPUTS:
- `energy_vector::Vector{Float32}` -> Energies of generated photons, in J
- `snames::NTuple{2, NTuple{3, String}}` -> Tuple congtaining 2 elements. 
    * Each of these elements is a 3D Tuple containing names of parent, and daughter (x2) species
    * The first Tuple will correspond to the ionisation, and will contain name of parent (e.g. H2O), ionized parent (H2O(+)), and electron (e(-)). NO ELECTRONIC STATES MUST BE GIVEN HERE IN PARENTHESIS
    * The second to the second consequent dissociation and will contain name of parent, heavy daughter species, and light daughter species. This will POTENTIALLY Includes ELECTRONIC STATES IN PARENTHESIS
- `tsh_energy::NTuple{2, Float32}` -> Ionisation energy for first reaction, and Dissociation Energy for total process, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 2-element vector. Each element is a NTuple{3, Float64} -> Each 3D Tuple contains the velocity vector for the 2 daughter species specified in product_types 
- `product_types` -> 2-element vector containing String objects -> Name of the 2 daughter species. 
    Only considers "plain" atomic/molecular species, meaning no information about electronic states is given.
"""
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
    call_photodestruction_logic_multiple(current_reaction, photochemical_info, parent_velocity, sun_tuple, photon_energy)

# OBJECTIVE:
- Call the apropriate photodestruction logic providing the relevant photochemical information

# INPUTS:
- `current_reaction::CurrentReaction` -> contains information on current photoreaction (see photodatabase.jl)
- `photochemical_info::Species` -> contains photochemical dtabase for parent species (see general_construct.jl)
- `parent_velocity::Float32` -> Parent velocity magnitude
- `sun_tuple:: NTuple{3, Float32}` -> 3D tuple containing direction of incident photons
- `energy_vector::Vector{Float32}` -> Energies of generated photons, in J
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