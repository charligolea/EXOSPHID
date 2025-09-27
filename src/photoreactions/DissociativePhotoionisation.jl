module DissociativePhotoionisation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

const c = 299792458        # Speed of light in m/s
const h  = 6.62607015e-34  # Plank Constant in J.s
const m_el = 9.1093837e-31 # Electron mass in kg
const m_fund = 1.66054e-27 # 1 M.U.


struct PhotoReaction
    E_bond_eV::Float64 # Total energy threshold for dissociative ionisation
    E_ionisation_eV::Float64 # Ionisation threshold energy in eV
    parent_velocity::Union{Tuple{Float64, Float64, Float64}, Float64} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::Union{Tuple{Float64, Float64, Float64}, Nothing}
    wvl_range::Union{Vector{Float64}, Nothing} # Contains photon wavelength for specified range. Not necessary if just running simple_photodissociation / simple_ionisation function
    energy_vector::Union{Vector{Float64}, Nothing} # Contains photon energies for specified range. Not necessary if just running simple_photodissociation / simple_ionisation function
    product_names::Union{Tuple{String, String, String}, Nothing} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::Union{Tuple{String, String, String}, Nothing}
    E_bond::Float64 # CALCULATED -> Total energy threshold for dissociative ionisation in J
    E_ionisation::Float64  # CALCULATED -> Ionisation energy in J
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_bond_eV, E_ionisation_eV, parent_velocity, sun_tuple, wvl_range, energy_vector, product_names, display_info)
        E_bond = E_bond_eV * 1.602e-19
        E_ionisation = E_ionisation_eV * 1.602e-19
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(E_bond_eV, E_ionisation_eV, parent_velocity, sun_tuple, wvl_range, energy_vector, product_names, product_types, E_bond, E_ionisation, display_info)
    end
end

function random_unit_tuple()
    θ, φ = 2π * rand(), π * rand()
    return (sin(φ) * cos(θ), sin(φ) * sin(θ), cos(φ))
end

function get_masses(parent_name, heavy_child_name, light_child_name)
    # Get masses for involved photoreaction
    # Nomenclature: Usually, water or hydrogen based photoreactions will result in a lighter product (like H, H2) and a heavier product (O, OH)

    possible_species = ("H", "H2", "O", "OH", "H2O", "HO2", "H2O2", "He", "Ne")
    mass_dict = (1* m_fund, 2 * m_fund, 16 * m_fund, 17 * m_fund, 18 * m_fund, 33 * m_fund, 34 * m_fund, 4 * m_fund, 20 * m_fund)

    m_parent = m_heavy = m_light = 0.0
    m_parent = mass_dict[findfirst(isequal(parent_name), possible_species)]
    m_heavy = mass_dict[findfirst(isequal(heavy_child_name), possible_species)]
    m_light = mass_dict[findfirst(isequal(light_child_name), possible_species)]

    return m_parent, m_heavy, m_light
end

function get_vibrorotational_energy(species::String)

    """
    Return the vibro-rotational energy [J] for the given `species`.

    Valid species values: "H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne".
    """
    energy = 0.0
    conversion_factor = 1.602e-19

    if species == "H2O"
        energy = 0.35 * conversion_factor
    elseif species == "O" || species == "O(3P)" || species == "O(1D)" || species == "O(1S)"
        energy = 0.0 * conversion_factor
    elseif species == "H2"
        energy = 0.30 * conversion_factor
    elseif species == "OH"
        energy = 0.20 * conversion_factor
    elseif species == "H(1s)" || species == "H(2s,2p)" || species == "H"
        energy = 0.0 * conversion_factor
    else
        throw(ArgumentError("Unknown species: $species"))
    end

    return energy
end

function allocate_velocity_ionisation(reaction, E_photon, p_photon)

    # Calculate velocities for the photoionisation products according to the moment and energy conservation quations
    # E_photon: Photon Energy in J
    # p_photon: linear momentum of photon in kg*m/s

    # 1. Get masses for all the species involved
    m_parent, m_heavy, m_light = get_masses(reaction.product_types[1], reaction.product_types[2], reaction.product_types[3])
    m_ion = m_parent

    # 2. Calculate excess energy for reaction
    E_excess = E_photon - reaction.E_ionisation

    # 3. Calculate unitary photon momentum vector
    u_ph = p_photon ./ norm(p_photon)

    # 4. Generate random direction for parent molecule and calculate velocity vector
    if reaction.parent_velocity isa Float64
        v_parent = reaction.parent_velocity .* random_unit_tuple()
    elseif reaction.parent_velocity isa Tuple{Float64, Float64, Float64}
        v_parent = reaction.parent_velocity
    end

    # 5. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ m_parent .* v_parent

    # 6. Calculate unitary vectors for electron and product ion
    u_el = random_unit_tuple()
    u_ion = random_unit_tuple()

    # 7. Calculate speed values for electron and product ion
    a = m_el * (m_ion * (dot(u_ion,u_ph)))^2 + (m_el^2 * m_ion)*(dot(u_el, u_ph)^2)
    b = -2 * (dot(u_el, u_ph) * m_el * m_ion) * (dot(total_momentum,u_ph))
    c = m_ion * dot(total_momentum,u_ph)^2 - 2 * E_excess * (m_ion * (dot(u_ion,u_ph)))^2

    Δ = b^2 - 4 * a * c
    v_el_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_el_tuple = v_el_mag .* u_el

    p_ion_tuple = (total_momentum .- m_el.*v_el_tuple)

    return p_ion_tuple
end

function allocate_velocity_dissociation(reaction, E_photon, p_photon, p_parent_ion)

    # Calculate velocities for the dissociation products according to the moment and energy conservation quations
    # E_photon: Photon Energy in J
    # p_photon: linear momentum of photon in kg*m/s
    # p_parent_ion: linear momentum of the previously ionised ion in kg m / s

    # 1. Get masses for all the species involved
    m_parent, m_heavy, m_light = get_masses(reaction.product_types[1], reaction.product_types[2], reaction.product_types[3])

    # 2. Calculate excess energy for reaction
    E_parent = 1/2 * m_parent * norm(reaction.parent_velocity)^2 + get_vibrorotational_energy(reaction.product_types[1])
    E_excess = (E_photon + E_parent) - (reaction.E_bond + get_vibrorotational_energy(reaction.product_types[2]) + get_vibrorotational_energy(reaction.product_types[3]))


    # 3. Calculate unitary parent ion vector
    u_parent_ion = p_parent_ion ./ norm(p_parent_ion)

    # 5. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ p_parent_ion

    # 6. Calculate unitary vectors for light and heavy products
    u_light = random_unit_tuple()
    u_heavy = random_unit_tuple()

    # 7. Calculate speed values for light and heavy products
    a = m_light * (m_heavy * (dot(u_heavy,u_parent_ion)))^2 + (m_light^2 * m_heavy)*(dot(u_light, u_parent_ion)^2)
    b = -2 * (dot(u_light, u_parent_ion) * m_light * m_heavy) * (dot(total_momentum,u_parent_ion))
    c = m_heavy * dot(total_momentum,u_parent_ion)^2 - 2 * E_excess * (m_heavy * (dot(u_heavy,u_parent_ion)))^2

    Δ = b^2 - 4 * a * c
    v_light_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_light_tuple = v_light_mag .* u_light

    v_heavy_tuple = (total_momentum .- m_light.*v_light_tuple) ./ m_heavy

    return v_light_tuple, v_heavy_tuple
end

# Simulate a single photodissociation reaction
function simulate_photoionisation(reaction::PhotoReaction, E_photon)

    # 1. Calculate photon linear momentum magnitude
    p_photon_magnitude = (E_photon) / c

    # 2. Calculate photon linear momentum vector by generating random incoming direction
    if reaction.sun_tuple isa Tuple{Float64, Float64, Float64}
        p_photon_unitary = sun_tuple
    elseif reaction.sun_tuple isa Nothing
        p_photon_unitary = random_unit_tuple()
    end

    p_photon = p_photon_magnitude .* p_photon_unitary

    # 3.1 Calculate linear momentum of parent ion after ionisation
    p_parent_ion = allocate_velocity_ionisation(reaction, E_photon, p_photon)

    # 3.2 Calculate velocity of heavy and light photoreaction product. NOTE that the photon has lost energy.
    p_photon_step2 = ((E_photon-reaction.E_ionisation) / c) .* p_photon_unitary
    v_light_tuple, v_heavy_tuple = allocate_velocity_dissociation(reaction, E_photon, p_photon_step2, p_parent_ion)

    return v_light_tuple, v_heavy_tuple
    
end

# Function to simulate multiple photodissociation reactions
function multiple_photoionisation(reaction::PhotoReaction)

    println("Simulating photoionisation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *""* " + " *""* reaction.product_names[3] *""* " + e-")

    final_speeds_light = []
    final_speeds_heavy = []

    # 1. Loop over photon energy vector
    for (index, E_photon) in enumerate(reaction.energy_vector)
        
        if E_photon > reaction.E_bond
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_light_tuple, v_heavy_tuple = simulate_photoionisation(reaction, E_photon)
            
            # 1.2. Store trajectories and speeds of every individual product
            push!(final_speeds_light, v_light_tuple./1000)
            push!(final_speeds_heavy, v_heavy_tuple./1000)

        end
    end

    final_speeds_heavy_norms = [norm(p) for p in final_speeds_heavy]
    final_speeds_light_norms = [norm(p) for p in final_speeds_light]


    # 4. Show mean, STD and median speeds for photo products
    if reaction.display_info
        data_speeds = DataFrame(
        Product = [reaction.product_names[2], reaction.product_names[3]],
        Mean_Speed = [mean(final_speeds_heavy_norms), mean(final_speeds_light_norms)],
        Median_Speed = [median(final_speeds_heavy_norms), median(final_speeds_light_norms)],
        STD_half = [std(final_speeds_heavy_norms)/2, std(final_speeds_light_norms)/2])
        println(data_speeds)
        println("")
    end

    return final_speeds_light, final_speeds_heavy
end

end

