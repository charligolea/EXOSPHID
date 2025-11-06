# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

"""
`exosphid_species`: atomic and moelcular species included in the EXOSPHID model
"""
const exosphid_species = ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne")

"""
`c`: Speed of light in m/s
"""
const c = 2.99792458f8               # Speed of light in m/s

"""
`m_el`: Electron mass in kg
"""
const m_el = 9.1093837e-31           # Electron mass in kg

"""
`m_fund`: 1 A.M.U.
"""
const m_fund = 1.66054e-27           # 1 A.M.U.

"""
`conversion_factor`: Conversion from eV to J
"""
const conversion_factor = 1.602f-19  # Conversion from eV to J

"""
`mass_species`: Relevant atomic / molecular species that appear either as parent or daughter products in the EXOSPHID database
"""
const mass_species = ("H", "H(-)", "H2", "O", "OH", "H2O", "HO2", "H2O2", "He", "Ne")

"""
`mass_dict`: Masses corresponding to the species in mass_species
"""
const mass_dict = (1* m_fund, 1* m_fund, 2 * m_fund, 16 * m_fund, 17 * m_fund, 18 * m_fund, 33 * m_fund, 34 * m_fund, 4 * m_fund, 20 * m_fund)


"""
CurrentReaction

- STRUCT Type object
- Stores the information for the photoreaction that is occuring for a specific iteration
- Includes:
    * `present_reaction::String` -> given reaction in Species.reaction_types
    * `reaction_name::String` -> given reaction in Species.reaction_names
    * `reaction_index::Int` -> Element of corresponding reaction in Species.reaction_types
    * `wvl_range::NTuple{2, Float32}` -> wavelength range in Species.reaction_types according to photon wavelength
"""
struct CurrentReaction
    present_reaction::String
    reaction_name::String
    reaction_index::Int
    wvl_range::NTuple{2, Float32}
end

# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

"""
    get_photodestruction_rates(species, solar_activity, dist_to_sun)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Return the photodestruction rate `k` [1/s] for a given `species` and `solar_activity` level (0 = quiet, 1 = active).

# INPUTS:
- `ph_info::Species` -> photochemical database of parent species
- `solar_activity::Float32` -> scalar from 0 (Quiet Sun) to 1 (Active Sun)
- `dist_to_sun::Float32` -> distance to Sun in A.U. (default for Moon = 1)
"""
function get_photodestruction_rates(ph_info::Species, solar_activity::Float32, dist_to_sun::Float32)

    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"
    @assert 0.0 <= dist_to_sun <= 100.0 "Distance to Sun must be positive and within known Solar System bounds (approx 100 A.U.)!"
    
    k = (1 - solar_activity) * ph_info.quiet_rate + solar_activity * ph_info.active_rate

    return k/(dist_to_sun^2)

end

get_photodestruction_rates(ph_info::Species, solar_activity::Real, dist_to_sun::Real) = get_photodestruction_rates(ph_info, Float32(solar_activity), Float32(dist_to_sun))


"""
    is_photoreaction_occuring(k, dt)
-------------------------------------------------------------------------------------------
# OBJECTIVE:
- Calculate probability of photoreaction from photodestruction rate
- If a generated random number between 0 and 1 is smaller than that probability, then the reaction will happening
- In that case, the function will return "true"

# INPUTS:
- `k::Float32` -> photodestruction rate
- `dt::Float32` -> time window
"""
function is_photoreaction_occuring(k::Float32, dt::Float32)
    @assert 0.0 <= k "Photodestruction rate must be positive!"
    @assert 0.0 <= dt "dt must be positive!"
    p = 1-exp(-k*dt)
    return rand() < p
end

is_photoreaction_occuring(k::Real, dt::Real) = is_photoreaction_occuring(Float32(k), Float32(dt))


"""
    get_wvl_threshold(ph_info)
-------------------------------------------------------------------------------------------
# OBJECTIVE:
- Get wavelength threshold in Angstrom for the photolysis of a certain species

# INPUTS:
- `ph_info::Species` -> Species type object containing photochemical info for given parent type

# REFERENCES:
- Huebner, W. F., Keady, J. J., & Lyon, S. P. (1992). Solar Photo Rates for Planetary Atmospheres and Atmospheric Pollutants. Astrophysics and Space Science, 195, 1–294.
- Huebner, W. F., & Mukherjee, J. (2015). Photoionization and photodissociation rates in solar and blackbody radiation fields. Planetary and Space Science, 106, 11–45.
"""
function get_wvl_threshold(ph_info::Species)
    return ph_info.wvl_threshold
end

"""
    get_species_photochemical_info(parent_type)
-------------------------------------------------------------------------------------------
# OBJECTIVE:
- Get photochemical data for the species to simulate

# INPUTS:
- `parent_type::String` -> name of parent species. Must be in `exosphid_species``

# OUTPUTS: Returns object of type Species containing the photochemical database for the given parent species, including:

- `tsh_energies`: Bond or ionization reaction. This will have the following types:
    * SPD: Single Float (dissociatione energy)
    * SPI: Single Float (ionization energy)
    * DPD: 2D Tuple (first dissociation energy, total dissociation energy)
    * DPI: 2D Tuple (first ionization energy, total ionization energy)
    * DiPI: 2D Tuple (ionization energy, dissociation energy)
- `species_names`:
    * SPD/SPI: 3D Tuple with parent, heavy product and light product species product_names
    * DPD/DPI/DiPI: 2D Tuple containing 2 3D Tuples for the sperate reactions
- `reaction_probabilities`: Branching ratio for a certain reaction in a certain wavelength range
- `reaction_names`: Identifiers for the different reactions
- `reaction_types`: 
    * SPD: Simple PhotoDissociation
    * SPI: Simple PhotoIonisation
    * DPD: Double PhotoDissociation 
    * DPI: Double PhotoIonisation
    * DiPI: Dissociative PhotoIonisation
- `wavelength_range`: wavelength ranges for which photochemical reactions are relevant for a given species

# REFERENCES:
- Combi, M. R., Harris, W. M., & Smyth, W. H. (2004). Gas Dynamics and Kinetics in the Cometary Coma: Theory and Observations.
- Huebner, W. F., Keady, J. J., & Lyon, S. P. (1992). Solar Photo Rates for Planetary Atmospheres and Atmospheric Pollutants. Astrophysics and Space Science, 195, 1–294.
- Huebner, W. F., & Mukherjee, J. (2015). Photoionization and photodissociation rates in solar and blackbody radiation fields. Planetary and Space Science, 106, 11–45.
- Marr, G.V. & West, J.B. (1976). Absolute Photoionization Cross-Section Tables for Helium, Neon, Argon and Krypton in the VUV Spectral Regions. Atomic Data and Nuclear Data Tables, 18, 497.
- Broad, John T. & Reinhardt, William P. (1976). One- and two-electron photoejection from H⁻: A multichannel J-matrix calculation. Physical Review A, 14(6), 2159–2173.
"""
function get_species_photochemical_info(parent_type::String)

    @assert parent_type in exosphid_species "Invalid parent species: $parent_type"

    if parent_type == "H2O"
        return H2O
    elseif parent_type == "OH"
        return OH
    elseif parent_type == "H2"
        return H2
    elseif parent_type == "H"
        return H
    elseif parent_type == "H(-)"
        return H_
    elseif parent_type == "HO2"
        return HO2
    elseif parent_type == "H2O2"
        return H2O2
    elseif parent_type == "He"
        return He
    elseif parent_type == "Ne"
        return Ne
    end
end


"""
    get_current_reaction(photon_wvl, photochemical_info)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Determine photoreaction that will occure and get photochemical info for such reaction

# INPUTS:
- `photochemical_info::Species` -> Contains photochemical database for parent species
- `photon_wvl::Real` -> scalar with photon wavelength in Angstrom to determine wavelength range

# OUTPUTS: CurrentReaction type object containing
- `present_reaction`: type of photodestruction process (SPD, SPI, DPD, DPI, DiPI, nothing)
- `reaction_name`: identifier
- `reaction_index`: reaction index within all the possible reactions for that aprent species
- `wvl_range`: in which wavelength range relevant for the parent species is the photon wavlength
end
"""
function get_current_reaction(photon_wvl::Real, photochemical_info::Species)

    wavelength_range = get_wavelength_range(photochemical_info)
    wvl_idx = findfirst(r -> r[1] < photon_wvl <= r[2], wavelength_range)

    if wvl_idx !== nothing
        wvl_range = wavelength_range[wvl_idx]
        r_idx = get_reaction_index(photochemical_info, wvl_idx)
        pr = get_reaction_type(photochemical_info, r_idx)
        rn = get_reaction_identifier(photochemical_info, r_idx)
        return CurrentReaction(pr, rn, r_idx, wvl_range)
    else
        return nothing, nothing, nothing, nothing
    end

end


"""
    get_reaction_index(ph_info, index)
    
# Objective:
- Subfunction of get_current_reaction
- for given wavelength range (identified by index) determine the reaction taht will happen according to the branchin ratios in that wavelength range

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
- `index::Int` -> Index corresponding to the current wavelength range of interest

"""
get_reaction_index(ph_info::Species, index::Int) = rand(Categorical([p[index] for p in ph_info.reaction_probabilities]))

"""
    get_reaction_identifier(photo_info, reaction_index)

# Objective:
- Subfunction of get_current_reaction
- get reaction identifier once the reaction that will happen is known from previous step (identified by reaction_index)

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
- `reaction_index::Int` -> Index corresponding to the current reaction
"""
get_reaction_identifier(ph_info::Species, reaction_index::Int) = ph_info.reaction_names[reaction_index]

"""
    get_reaction_type(photo_info, reaction_index): 

# Objective:
- Subfunction of get_current_reaction
- find reaction type (SPD, SPI...etc) once the reaction that will happen is known from previous step (identified by reaction_index)

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
- `reaction_index::Int` -> Index corresponding to the current reaction
"""
get_reaction_type(ph_info::Species, reaction_index::Int) = ph_info.reaction_types[reaction_index]

"""
    get_wavelength_range(photo_info)

# Objective:
- Subfunction of get_current_reaction
- unpack wavelength ranges for given species

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
"""
get_wavelength_range(ph_info::Species) = ph_info.wavelength_range


"""
    get_vibrorotational_energy(species)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Return the vibro-rotational energy [J] for a given species

# INPUTS:
- `species::String` -> Species identifier including electronic state
    Valid values: "H2O", "OH", "OH(X2π)", "OH(A2Σ+)", "OH(1Σ+)", "OH(12Δ/22Π)", 
                  "OH(A2Σ+, v'=2)", "OH(A2Σ+, v'=3)", "OH(B2Σ)", "OH(D2Σ)", 
                  "O", "O(3P)", "O(1D)", "O(1S)", "H2", "H(1s)", "H(2s,2p)", "H", 
                  "H(-)", "HO2", "H2O2", "He", "Ne"

# OUTPUTS:
- `energy::Float32` -> Vibro-rotational energy in Joules
"""
function get_vibrorotational_energy(species::String)

    if species == "H2O"
        energy = 0.35f0 * conversion_factor
    elseif species == "OH(X2π)" || species == "OH"
        energy = 0.22f0 * conversion_factor
    elseif species == "OH(A2Σ+)"
        energy = 0.25f0 * conversion_factor
    elseif species == "OH(DPD)"
        energy = 0.00f0 * conversion_factor
    elseif species == "OH(1Σ+)" || species == "OH(12Δ/22Π)"
        energy = 0.25f0 * conversion_factor
    elseif species == "OH(A2Σ+, v'=3)"
        # energy = 1.50f0 * conversion_factor
        energy = 0.20f0 * conversion_factor # Considered as ground state, see vibrorotational theory
    elseif species == "OH(A2Σ+, v'=2)"
        energy = 1.05f0 * conversion_factor
    elseif species == "OH(B2Σ)"
        energy = 0.05f0 * conversion_factor
    elseif species == "OH(D2Σ)"
        energy = 0.20f0 * conversion_factor
    elseif species == "H2"
        energy = 0.30f0 * conversion_factor
    elseif species in ("O", "O(3P)", "O(1D)", "O(1S)", "H(1s)", "H(2s,2p)" , "H" , "H(-)" , "He", "Ne")
        energy = 0.0f0 * conversion_factor
    elseif species == "HO2"
        energy = 0.43f0 * conversion_factor
    elseif species == "H2O2"
        energy = 0.45f0 * conversion_factor
    else
        energy = nothing
    end

    return energy
end


"""
    get_electronic_energy_predis(species)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Return the electronic energy [J] for a species in predissociation cases (OH)

# INPUTS:
- `species::String` -> Species identifier with electronic state
    Valid values: "OH(X2π)", "OH", "OH(A2Σ+)", "OH(A2Σ+, v'=2)", "OH(A2Σ+, v'=3)", 
                 "OH(1Σ+)", "OH(12Δ/22Π)", "OH(B2Σ)", "OH(D2Σ)"

# OUTPUTS:
- `energy::Float32` : Electronic energy in Joules
"""
function get_electronic_energy_predis(species::String)

    if species == "OH(X2π)" || species == "OH"
        energy = 0.00f0 * conversion_factor
    elseif species == "OH(A2Σ+)" || species == "OH(A2Σ+, v'=3)" || species == "OH(A2Σ+, v'=2)"
        energy = 4.05f0 * conversion_factor
    elseif species == "OH(1Σ+)"
        energy = 4.05f0 * conversion_factor
    elseif species == "OH(12Δ/22Π)"
        energy = 6.50f0 * conversion_factor
    elseif species == "OH(B2Σ)"
        energy = 8.65f0 * conversion_factor
    elseif species == "OH(D2Σ)"
        energy = 10.18f0 * conversion_factor
    else
        energy = nothing
    end

    return energy
end


"""
    get_masses(parent_name; heavy_child_name=nothing, light_child_name=nothing, mode="PD")
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Get masses for the parent and child species involved in a photoreaction
- For PD (photodissociation): returns masses of parent, heavy child, and light child
- For PI (photoionization): returns mass of parent only

# INPUTS:
- `parent_name::String` -> Name of parent species
- `heavy_child_name::Union{String, Nothing}` -> Name of heavier child species (required for PD)
- `light_child_name::Union{String, Nothing}` -> Name of lighter child species (required for PD)
- `mode::String` -> "PD" for photodissociation, "PI" for photoionization
"""
function get_masses(parent_name; heavy_child_name=nothing, light_child_name=nothing, mode="PD")

    m_parent = mass_dict[findfirst(isequal(parent_name), mass_species)]

    if mode == "PI"
        return m_parent
    elseif mode == "PD"
        @assert heavy_child_name !== nothing "heavy_child_name must be provided for PD"
        @assert light_child_name !== nothing "light_child_name must be provided for PD"
        m_heavy = mass_dict[findfirst(isequal(heavy_child_name), mass_species)]
        m_light = mass_dict[findfirst(isequal(light_child_name), mass_species)]
        return m_parent, m_heavy, m_light
    end
end


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

export get_species_photochemical_info
export get_current_reaction
export exosphid_species