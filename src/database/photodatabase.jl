# ─────────────────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────────────────

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
`mass_species`: Relevant atomic / molecular species that appear either as parent or daughter 
    products in the EXOSPHID database
"""
const mass_species = ("H", "H(-)", "H2", "O", "OH", "H2O", "HO2", "H2O2", "He", "Ne")


"""
CurrentReaction

- STRUCT Type object
- Stores the information for the photoreaction that is occuring for a specific iteration
- Includes:
    * `present_reaction::String` -> given reaction in Species.reaction_types
    * `reaction_name::String` -> given reaction in Species.reaction_names
    * `reaction_index::Integer` -> Element of corresponding reaction in 
        Species.reaction_types
    * `wvl_range::NTuple{2, Float32}` -> wavelength range in Species.reaction_types 
        according to photon wavelength
"""
struct CurrentReaction
    present_reaction::String
    reaction_name::String
    reaction_index::Integer
    wvl_range::NTuple{2, Float32}
end

# ─────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────────────

"""
eV2J(eV)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Conversion from eV to J

# INPUTS:
- `eV::Real` -> Energy in eV
"""
eV2J(eV::Float32) = 1.602f-19 * eV
eV2J(eV::Real) = eV2J(Float32(eV))

"""
amu2kg(m)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Conversion from amu (atomic mass units) to kg

# INPUTS:
- `m::Real` -> Mass in amu
"""
amu2kg(m::Integer) = 1.66054e-27 * m
amu2kg(m::Real) = amu2kg(Int(m))

"""
`mass_dict`: Masses corresponding to the species in mass_species
"""
const mass_dict = (amu2kg(1), amu2kg(1), amu2kg(2), amu2kg(16), amu2kg(17), 
                   amu2kg(18), amu2kg(33), amu2kg(34), amu2kg(4), amu2kg(20))

"""
    get_photodestruction_rates(species, solar_activity, dist_to_sun)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Return the photodestruction rate `k` [1/s] for a given `species` and `solar_activity` 
    level (0 = quiet, 1 = active).

# INPUTS:
- `ph_info::Species` -> photochemical database of parent species
- `solar_activity::Float32` -> scalar from 0 (Quiet Sun) to 1 (Active Sun)
- `dist_to_sun::Float32` -> distance to Sun in A.U. (default for Moon = 1)
"""
function get_photodestruction_rates(ph_info::Species, solar_activity::Float32, 
dist_to_sun::Float32)

    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"
    @assert 0.0 <= dist_to_sun <= 100.0 "Distance to Sun must be positive and within known 
        Solar System bounds (approx 100 A.U.)!"
    
    k = (1 - solar_activity) * ph_info.quiet_rate + solar_activity * ph_info.active_rate

    return k/(dist_to_sun^2)

end

function get_photodestruction_rates(ph_info::Species, solar_activity::Real, dist_to_sun::Real)
    return get_photodestruction_rates(ph_info, Float32(solar_activity), Float32(dist_to_sun))
end


"""
    is_photoreaction_occuring(k, dt)
-------------------------------------------------------------------------------------------
# OBJECTIVE:
- Calculate probability of photoreaction from photodestruction rate
- If a generated random number between 0 and 1 is smaller than that probability, then the 
    reaction will happening
- In that case, the function will return "true"

# INPUTS:
- `k::Float32` -> photodestruction rate
- `dt::Float32` -> time window
"""
function is_photoreaction_occuring(k::Real, dt::Real)
    @assert 0.0 <= k "Photodestruction rate must be positive!"
    @assert 0.0 <= dt "dt must be positive!"
    p = 1-exp(-k*dt)
    return rand() < p
end


"""
    get_wvl_threshold(ph_info)
-------------------------------------------------------------------------------------------
# OBJECTIVE:
- Get wavelength threshold in Angstrom for the photolysis of a certain species

# INPUTS:
- `ph_info::Species` -> Species type object containing photochemical info for given parent 
    type

# REFERENCES:
- Huebner, W. F., Keady, J. J., & Lyon, S. P. (1992). Solar Photo Rates for Planetary 
    Atmospheres and Atmospheric Pollutants. Astrophysics and Space Science, 195, 1–294.
- Huebner, W. F., & Mukherjee, J. (2015). Photoionization and photodissociation rates in 
    solar and blackbody radiation fields. Planetary and Space Science, 106, 11–45.
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

# OUTPUTS: Returns object of type Species containing the photochemical database for the 
    given parent species, including:

- `tsh_energies`: Bond or ionization reaction. This will have the following types:
    * SPD: Single Float (dissociatione energy)
    * SPI: Single Float (ionization energy)
    * DPD: 2D Tuple (first dissociation energy, total dissociation energy)
    * DPI: 2D Tuple (first ionization energy, total ionization energy)
    * DiPI: 2D Tuple (ionization energy, dissociation energy)
- `species_names`:
    * SPD/SPI: 3D Tuple with parent, heavy product and light product species product_names
    * DPD/DPI/DiPI: 2D Tuple containing 2 3D Tuples for the sperate reactions
- `reaction_probabilities`: Branching ratio for a certain reaction in a certain wavelength 
    range
- `reaction_names`: Identifiers for the different reactions
- `reaction_types`: 
    * SPD: Simple PhotoDissociation
    * SPI: Simple PhotoIonisation
    * DPD: Double PhotoDissociation 
    * DPI: Double PhotoIonisation
    * DiPI: Dissociative PhotoIonisation
- `wavelength_range`: wavelength ranges for which photochemical reactions are relevant for a 
    given species

# REFERENCES:
- Combi, M. R., Harris, W. M., & Smyth, W. H. (2004). Gas Dynamics and Kinetics in the 
    Cometary Coma: Theory and Observations.
- Huebner, W. F., Keady, J. J., & Lyon, S. P. (1992). Solar Photo Rates for Planetary 
    Atmospheres and Atmospheric Pollutants. Astrophysics and Space Science, 195, 1–294.
- Huebner, W. F., & Mukherjee, J. (2015). Photoionization and photodissociation rates in 
    solar and blackbody radiation fields. Planetary and Space Science, 106, 11–45.
- Marr, G.V. & West, J.B. (1976). Absolute Photoionization Cross-Section Tables for Helium, 
    Neon, Argon and Krypton in the VUV Spectral Regions. Atomic Data and Nuclear Data Tables, 18, 497.
- Broad, John T. & Reinhardt, William P. (1976). One- and two-electron photoejection from H⁻: 
    A multichannel J-matrix calculation. Physical Review A, 14(6), 2159–2173.
"""
function get_species_photochemical_info(parent_type::String)

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
    throw(ArgumentError("Invalid parent species: $parent_type"))
end


"""
    get_current_reaction(photon_wvl, photochemical_info)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Determine photoreaction that will occure and get photochemical info for such reaction

# INPUTS:
- `photochemical_info::Species` -> Contains photochemical database for parent species
- `photon_wvl::Real` -> scalar with photon wavelength in Angstrom to determine wavelength 
    range

# OUTPUTS: CurrentReaction type object containing
- `present_reaction`: type of photodestruction process (SPD, SPI, DPD, DPI, DiPI, nothing)
- `reaction_name`: identifier
- `reaction_index`: reaction index within all the possible reactions for that aprent species
- `wvl_range`: in which wavelength range relevant for the parent species is the photon 
    wavlength
"""
function get_current_reaction(photon_wvl::Real, photochemical_info::Species)

    wavelength_range = get_wavelength_range(photochemical_info)
    wvl_idx = findfirst(r -> r[1] < photon_wvl <= r[2], wavelength_range)

    if !isnothing(wvl_idx)
        wvl_range = wavelength_range[wvl_idx]
        r_idx = get_reaction_index(photochemical_info, wvl_idx)
        pr = get_reaction_type(photochemical_info, r_idx)
        rn = get_reaction_identifier(photochemical_info, r_idx)
        return CurrentReaction(pr, rn, r_idx, wvl_range)
    else
        return "", "", Integer[], ()
    end

end


"""
    get_reaction_index(ph_info, index)
-------------------------------------------------------------------------------------------
    
# Objective:
- Subfunction of get_current_reaction
- for given wavelength range (identified by index) determine the reaction taht will happen 
    according to the branchin ratios in that wavelength range

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
- `index::Integer` -> Index corresponding to the current wavelength range of interest

"""
function get_reaction_index(ph_info::Species, index::Integer)
    return rand(Categorical([p[index] for p in ph_info.reaction_probabilities]))
end


"""
    get_reaction_identifier(photo_info, reaction_index)
-------------------------------------------------------------------------------------------

# Objective:
- Subfunction of get_current_reaction
- get reaction identifier once the reaction that will happen is known from previous step 
    (identified by reaction_index)

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
- `reaction_index::Integer` -> Index corresponding to the current reaction
"""
function get_reaction_identifier(ph_info::Species, reaction_index::Integer)
    return ph_info.reaction_names[reaction_index]
end


"""
    get_reaction_type(photo_info, reaction_index): 
-------------------------------------------------------------------------------------------

# Objective:
- Subfunction of get_current_reaction
- find reaction type (SPD, SPI...etc) once the reaction that will happen is known from 
    previous step (identified by reaction_index)

# Inputs:
- `ph_info::Species` -> Contains photochemical database for parent species
- `reaction_index::Integer` -> Index corresponding to the current reaction
"""
function get_reaction_type(ph_info::Species, reaction_index::Integer)
    return ph_info.reaction_types[reaction_index]
end


"""
    get_wavelength_range(photo_info)
-------------------------------------------------------------------------------------------

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

    energy = 0.00f0

    if species == "H2O"
        energy = eV2J(0.35f0)
    elseif species == "OH(X2π)" || species == "OH"
        energy =eV2J(0.22f0)
    elseif species == "OH(A2Σ+)"
        energy = eV2J(0.25f0)
    elseif species == "OH(DPD)"
        energy = eV2J(0.00f0)
    elseif species == "OH(1Σ+)" || species == "OH(12Δ/22Π)"
        energy = eV2J(0.25f0)
    elseif species == "OH(A2Σ+, v'=3)"
        # energy = eV2J(1.50f0)
        energy = eV2J(0.20f0) # Considered as ground state, see vibrorotational theory
    elseif species == "OH(A2Σ+, v'=2)"
        energy =eV2J(1.05f0)
    elseif species == "OH(B2Σ)"
        energy = eV2J(0.05f0)
    elseif species == "OH(D2Σ)"
        energy = eV2J(0.20f0)
    elseif species == "H2"
        energy = eV2J(0.30f0)
    elseif species in ("O", "O(3P)", "O(1D)", "O(1S)", 
                       "H(1s)", "H(2s,2p)" , "H" , "H(-)",
                       "He", "Ne")
        energy = eV2J(0.0f0)
    elseif species == "HO2"
        energy = eV2J(0.43f0)
    elseif species == "H2O2"
        energy = eV2J(0.45f0)
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

    energy = 0.00f0

    if species == "OH(X2π)" || species == "OH"
        energy = eV2J(0.00f0)
    elseif species == "OH(A2Σ+)" || species == "OH(A2Σ+, v'=3)" || species == "OH(A2Σ+, v'=2)"
        energy = eV2J(4.05f0)
    elseif species == "OH(1Σ+)"
        energy = eV2J(4.05f0)
    elseif species == "OH(12Δ/22Π)"
        energy = eV2J(6.50f0)
    elseif species == "OH(B2Σ)"
        energy = eV2J(8.65f0)
    elseif species == "OH(D2Σ)"
        energy = eV2J(10.18f0)        
    end

    return energy
end


"""
    get_masses(species_name)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Get masses for a certain species

# INPUTS:
- `species_name`: nmame of species in question, without electronic state or charge information

# EXAMPLE:
- `H2O`, 'OH' are allowed
- `H2O(+)`, 'O(3P)' are `NOT` allowed
"""
get_masses(species_name) = mass_dict[findfirst(isequal(species_name), mass_species)]


# ─────────────────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ─────────────────────────────────────────────────────────────────────────────────────────

export get_species_photochemical_info
export get_current_reaction
export exosphid_species