# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS TO RUN PHOTOLYSIS LOGIC DEPEND ON REACTION TYPE
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

"""
    SPD_logic(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Handles Simple PhotoDissociation (SPD) logic

# INPUTS:
- `photon_energy::Float32` -> Photon Energy in J
- `snames::NTuple{3, String}` -> 3D Tuple containing names of parent, and daughter (x2) 
    species, POTENTIALLY INCLUDING ELECTRONIC STATES IN PARENTHESIS
- `tsh_energy::Float32` -> Bond/Dissociation Energy for parent species, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 2-element vector. Each element is a NTuple{3, Float64} -> Each 3D 
    Tuple contains the velocity vector for the 2 daughter species specified in product_types 
- `product_types` -> 2-element vector containing String objects -> Name of the 2 daughter 
    species. 
    Only considers "plain" atomic/molecular species, meaning no information about electronic 
        states is given.
"""
function SPD_logic(photon_energy::Float32, snames::NTuple{3, String}, tsh_energy::Float32, 
                parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), snames)
    reaction = PhotoReaction(eV2J(tsh_energy), parent_velocity, sun_tuple, snames, false)
    final_speeds_light, final_speeds_heavy = simulate_photodissociation(reaction, photon_energy)

    product_types = [product_names[2], product_names[3]]
    product_velocities = [final_speeds_heavy, final_speeds_light]

    return product_velocities, product_types
end

function SPD_logic(photon_energy::Real, snames::NTuple{3, String}, tsh_energy::Float32, 
                parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    return SPD_logic(Float32(photon_energy), snames, tsh_energy, parent_velocity, sun_tuple)
end


"""
    DPD_logic(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Handles Double PhotoDissociation (DPD) logic

# INPUTS:
- `photon_energy::Float32` -> Photon Energy in J
- `snames::NTuple{2, NTuple{3, String}}` -> Tuple congtaining 2 elements. 
    * Each of these elements is a 3D Tuple containing names of parent, and daughter (x2) 
        species
    * The first Tuple will corrspond to the first dissociation, the other to the second 
        consequent dissociation
    * POTENTIALLY INCLUDING ELECTRONIC STATES IN PARENTHESIS
- `tsh_energy::NTuple{2, Float32}` -> Bond/Dissociation Energy for first and for the total 
    reaction, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 3-element vector. Each element is a NTuple{3, Float64} -> Each 3D 
    Tuple contains the velocity vector for the 3 daughter species specified in product_types 
- `product_types` -> 3-element vector containing String objects -> Name of the 3 daughter 
    species. 
    * Only considers "plain" atomic/molecular species, meaning no information about electronic 
        states is given.
    * The first species will be the heavy product, the other 2 the lighter ones.
    * Currently, this reaction is only supported for the DPD of Water: 
        H2O -> OH + H -> O + H + H
"""
function DPD_logic(photon_energy::Float32, snames::NTuple{2, NTuple{3, String}}, 
                tsh_energy::NTuple{2, Float32}, parent_velocity::NTuple{3, Float32}, 
                sun_tuple::NTuple{3, Float32})
    
    product_names = (map(s -> replace(s, r"\((?!\+).*?\)" => ""), snames[1]), 
                    map(s -> replace(s, r"\((?!\+).*?\)" => ""), snames[2]))

    reaction = PhotoReaction(eV2J(tsh_energy[1]), parent_velocity, sun_tuple, snames[1], false)
    final_speeds_light_1, final_speeds_heavy_old = 
        simulate_photodissociation(reaction, photon_energy)

    reaction = PhotoReaction(eV2J(tsh_energy[2]), map(Float32, final_speeds_heavy_old), 
                            sun_tuple, snames[2], false)
    final_speeds_light_2, final_speeds_heavy = 
        simulate_photodissociation(reaction, photon_energy)
    
    product_types = [product_names[2][2], product_names[1][3], product_names[2][3]]
    product_velocities = [final_speeds_heavy, final_speeds_light_1, final_speeds_light_2]

    return product_velocities, product_types
end

function DPD_logic(photon_energy::Real, snames::NTuple{3, String}, tsh_energy::Float32, 
                parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    return DPD_logic(Float32(photon_energy), snames, tsh_energy, parent_velocity, sun_tuple)
end


"""
    SPI_logic(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Handles Simple PhotoIonisation (SPI) logic

# INPUTS:
- `photon_energy::Float32` -> Photon Energy in J
- `snames::NTuple{3, String}` -> 3D Tuple containing names of parent (e.g H2O), ionised 
    product (e.g. H2O(+)) and electron (e.g. e(-))
    * electronic states should not be included in this case
- `tsh_energy::Float32` -> Ionisation Energy in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 1-element vector, containing NTuple{3, Float64} -> The 3D Tuple 
    contains the velocity vector for the ionised parent species
- `product_types` -> 1-element vector containing String object -> Name of the ionised parent
"""
function SPI_logic(photon_energy::Float32, snames::NTuple{3, String}, tsh_energy::Float32, 
                parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    reaction = PhotoReaction(eV2J(tsh_energy), parent_velocity, sun_tuple, snames, false)
    product_velocities = [simulate_photoionisation(reaction, photon_energy)]
    product_types = [snames[2]]

    return product_velocities, product_types
end

function SPI_logic(photon_energy::Real, snames::NTuple{3, String}, tsh_energy::Float32, 
                parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    return SPI_logic(Float32(photon_energy), snames, tsh_energy, parent_velocity, sun_tuple)
end


"""
    DPI_logic(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Handles Double PhotoIonisation (DPI) logic

# INPUTS:
- `photon_energy::Float32` -> Photon Energy in J
- `snames::NTuple{3, String}` -> Tuple containing 2 Tuple elements.
    * Both of these elements are 3D Tuples containing names of parent (e.g H2O), ionised 
        product (e.g. H2O(+)) and electron (e.g. e(-))
    * Two tuples, for 2 ionisation reactions
    * electronic states should not be included in parenthesis in this case
- `tsh_energy::NTuple{2, Float32}` -> Ionisation Energy for first and for total reaction, 
    in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 1-element vector, containing NTuple{3, Float64} -> The 3D Tuple 
    contains the velocity vector for the final ionised product
- `product_types` -> 1-element vector containing String object -> Name of the final ionised 
    species
"""
function DPI_logic(photon_energy::Float32, snames::NTuple{2, NTuple{3, String}}, 
                tsh_energy::NTuple{2, Float32}, parent_velocity::NTuple{3, Float32}, 
                sun_tuple::NTuple{3, Float32})
    reaction = 
        PhotoReaction(eV2J(tsh_energy[1]), parent_velocity, sun_tuple, snames[1], false)
    final_speeds_ion = simulate_photoionisation(reaction, photon_energy)

    reaction = PhotoReaction(eV2J(tsh_energy[2]), map(Float32, final_speeds_ion), sun_tuple, 
                            snames[2], false)
    product_velocities = [simulate_photoionisation(reaction, photon_energy)]
    product_types = [snames[2][2]]

    return product_velocities, product_types
end

function DPI_logic(photon_energy::Real, snames::NTuple{3, String}, tsh_energy::Float32, 
                parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    return DPI_logic(Float32(photon_energy), snames, tsh_energy, parent_velocity, sun_tuple)
end


"""
    DiPI_logic(photon_energy, snames, tsh_energy, parent_velocity, sun_tuple): 
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Handles Dissociative PhotoIonisation (DiPI) logic
- The parent species is ionised, then dissociated.

# INPUTS:
- `photon_energy::Float32` -> Photon Energy in J
- `snames::NTuple{2, NTuple{3, String}}` -> Tuple congtaining 2 elements. 
    * Each of these elements is a 3D Tuple containing names of parent, and daughter (x2) 
        species
    * The first Tuple will correspond to the ionisation, and will contain name of parent 
        (e.g. H2O), ionized parent (H2O(+)), and electron (e(-)). 
        NO ELECTRONIC STATES MUST BE GIVEN HERE IN PARENTHESIS
    * The second to the second consequent dissociation and will contain name of parent, 
        heavy daughter species, and light daughter species. 
        This will POTENTIALLY Includes ELECTRONIC STATES IN PARENTHESIS
- `tsh_energy::NTuple{2, Float32}` -> Ionisation energy for first reaction, and 
    Dissociation Energy for total process, in J
- `parent_velocity::NTuple{3, Float32}` -> 3D velocity Tuple
- `sun_tuple::NTuple{3, Float32})` -> 3D Solar Tuple, or direction of incoming photon

# OUTPUTS:
- `product_velocities` -> 2-element vector. Each element is a NTuple{3, Float64} -> Each 3D 
    Tuple contains the velocity vector for the 2 daughter species specified in product_types 
- `product_types` -> 2-element vector containing String objects -> Name of the 2 daughter 
    species. Only considers "plain" atomic/molecular species, meaning no information about 
    electronic states is given.
"""
function DiPI_logic(photon_energy::Float32, snames::NTuple{2, NTuple{3, String}}, 
                tsh_energy::NTuple{2, Float32}, parent_velocity::NTuple{3, Float32}, 
                sun_tuple::NTuple{3, Float32})
    reaction = 
        PhotoReaction(eV2J(tsh_energy[1]), parent_velocity, sun_tuple, snames[1], false)
    final_speeds_ion = simulate_photoionisation(reaction, photon_energy)

    product_names = map(s -> replace(s, r"\(\+\)" => ""), snames[2])
    
    reaction = PhotoReaction(eV2J(tsh_energy[2]), map(Float32, final_speeds_ion), 
                            sun_tuple, product_names, false)
    final_speeds_light, final_speeds_heavy = 
        simulate_photodissociation(reaction, photon_energy)
    
    product_velocities = [final_speeds_heavy, final_speeds_light]
    product_types = [snames[2][2], snames[2][3]]

    return product_velocities, product_types
end

function DiPI_logic(photon_energy::Real, snames::NTuple{3, String}, tsh_energy::Float32, 
                    parent_velocity::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32})
    return DiPI_logic(Float32(photon_energy), snames, tsh_energy, 
        parent_velocity, sun_tuple)
end


"""
    call_photodestruction_logic(current_reaction, photochemical_info, parent_velocity, 
                                sun_tuple, photon_energy)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Call the apropriate photodestruction logic providing the relevant photochemical 
    information

# INPUTS:
- `current_reaction::CurrentReaction` -> contains information on current photoreaction 
    (see `photodatabase.jl`)
- `photochemical_info::Species` -> contains photochemical dtabase for parent species 
    (see `general_construct.jl`)
- `parent_velocity::NTuple{3, Float32}` -> 3D tuple with parent velocity
- `sun_tuple:: NTuple{3, Float32}` -> 3D tuple containing direction of incident photons
- `photon_energy::Real` -> in J
"""
function call_photodestruction_logic(current_reaction::CurrentReaction, 
                                    photochemical_info::Species, 
                                    parent_velocity::NTuple{3, Float32}, 
                                    sun_tuple::NTuple{3, Float32}, 
                                    photon_energy::Real)

    pr = current_reaction.present_reaction
    r_idx = current_reaction.reaction_index

    tsh_energies = get_tsh_energies(photochemical_info)
    species_names = get_species_names(photochemical_info)

    # Simple photodissociations
    if pr == "SPD"
        return SPD_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], 
                        parent_velocity, sun_tuple)
    
    # Double photodissociation
    elseif pr == "DPD" 
        return DPD_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], 
                        parent_velocity, sun_tuple)
    
    # Simple photoionisations
    elseif pr == "SPI"
        return SPI_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], 
                        parent_velocity, sun_tuple)

    # Double electron ejection of negative atomic hydrogen
    elseif pr == "DPI"
        return DPI_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], 
                        parent_velocity, sun_tuple)

    # Dissociative photoionisations
    elseif pr == "DiPI"
        return DiPI_logic(photon_energy, species_names[r_idx], tsh_energies[r_idx], 
                        parent_velocity, sun_tuple)
    end

    throw(ArgumentError("Invalid reaction type: $pr. 
                        Reaction type should be either SPD, SPI, DPD, DPI, DiPI"))

end


# ─────────────────────────────────────────────────────────────────────────────────────────
# SUBFUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────────────

"""
    get_tsh_energies(photo_info) 
-------------------------------------------------------------------------------------------

# OBJECTIVE: 
- unpack dissociation / ionisation energies for species of interest
"""
get_tsh_energies(photo_info) = photo_info.tsh_energies  # eV


"""
    get_species_names(photo_info) 
-------------------------------------------------------------------------------------------

# OBJECTIVE: 
- unpack parent and product species names for reactions for species of interest
"""
get_species_names(photo_info) = photo_info.species_names


# ─────────────────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ─────────────────────────────────────────────────────────────────────────────────────────
export call_photodestruction_logic