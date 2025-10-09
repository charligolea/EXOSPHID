module SimplePhotodissociation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

const c = 2.99792458f8        # Speed of light in m/s
const m_fund = 1.66054e-27 # 1 M.U.

struct PhotoReaction
    E_bond::Float32 # Threshold energy for given photodissociation reaction in J -> USER INPUT
    v_parent::NTuple{3, Float32} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::NTuple{3, Float32} # 
    product_names::NTuple{3, String} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::NTuple{3, String}
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_bond::Real, v_parent::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_bond), Float32.(v_parent), Float32.(sun_tuple), product_names, product_types, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::Real, sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        vp = Float32(v_parent) .* random_unit_tuple()
        new(Float32(E_bond), vp, Float32.(sun_tuple), product_names, product_types, display_info)
    end
    
    function PhotoReaction(E_bond::Real, v_parent::NTuple{3, Float32}, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_bond), Float32.(v_parent), random_unit_tuple(), product_names, product_types, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::Real, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        vp = Float32(v_parent) .* random_unit_tuple()
        new(Float32(E_bond), vp, random_unit_tuple(), product_names, product_types, display_info)
    end

end

function random_unit_tuple()
    θ, φ = Float32(2π) * rand(Float32), Float32(π) * rand(Float32)
    return (Float32(sin(φ) * cos(θ)), Float32(sin(φ) * sin(θ)), Float32(cos(φ)))
end

function get_masses(parent_name, heavy_child_name, light_child_name)
    # Get masses for involved photoreaction
    # Nomenclature: Usually, water or hydrogen based photoreactions will result in a lighter product (like H, H2) and a heavier product (O, OH)

    possible_species = ("H", "H2", "O", "OH", "H2O", "HO2", "H2O2", "He", "Ne")
    mass_dict = (1* m_fund, 2 * m_fund, 16 * m_fund, 17 * m_fund, 18 * m_fund, 33 * m_fund, 34 * m_fund, 4 * m_fund, 20 * m_fund)

    m_parent = m_heavy = m_light = 0.0f0
    m_parent = mass_dict[findfirst(isequal(parent_name), possible_species)]
    m_heavy = mass_dict[findfirst(isequal(heavy_child_name), possible_species)]
    m_light = mass_dict[findfirst(isequal(light_child_name), possible_species)]

    return m_parent, m_heavy, m_light
end

calculate_photon_momentum(E_photon, sun_tuple) = (E_photon/c) .* sun_tuple

function calculate_excess_energy(reaction::PhotoReaction, m_parent, E_photon)

    if reaction.product_types[1] == "OH" && reaction.product_names[1] != "OH(DPD)"
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2
        E_excess = (E_parent + E_photon - get_electronic_energy_predis(reaction.product_names[1]) - get_vibrorotational_energy(reaction.product_names[1]))
    elseif reaction.product_names[1] == "OH(DPD)"
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2
        E_excess = (E_photon + E_parent) - reaction.E_bond 
    else
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2 + get_vibrorotational_energy(reaction.product_names[1])
        E_excess = (E_photon + E_parent) - reaction.E_bond - (get_vibrorotational_energy(reaction.product_names[2]) + get_vibrorotational_energy(reaction.product_names[3]))
    end

    """   if E_excess <0
        print("E_excess: ", E_excess/1.602f-19, " eV\n")
    end"""

    return E_excess
end

function get_vibrorotational_energy(species::String)

    """
    Return the vibro-rotational energy [J] for the given `species`.

    Valid species values: "H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne".
    """
    energy = 0.0f0
    conversion_factor = 1.602f-19

    if species == "H2O"
        energy = 0.35f0 * conversion_factor
    elseif species == "OH(X2π)" || species == "OH"
        energy = 0.20f0 * conversion_factor
    elseif species == "OH(A2Σ+)"
        energy = 0.25f0 * conversion_factor
    elseif species == "OH(1Σ+)" || species == "OH(12Δ/22Π)"
        energy = 0.25f0 * conversion_factor
    elseif species == "OH(A2Σ+, v'=3)"
        # energy = 1.50f0 * conversion_factor
        energy = 0.20f0 * conversion_factor # TO REVIEW
    elseif species == "OH(A2Σ+, v'=2)"
        energy = 1.05f0 * conversion_factor
    elseif species == "OH(B2Σ)"
        energy = 0.05f0 * conversion_factor
    elseif species == "OH(D2Σ)"
        energy = 0.20f0 * conversion_factor
    elseif species == "O" || species == "O(3P)" || species == "O(1D)" || species == "O(1S)"
        energy = 0.0f0 * conversion_factor
    elseif species == "H2"
        energy = 0.30f0 * conversion_factor
    elseif species == "H(1s)" || species == "H(2s,2p)" || species == "H" || species == "H(-)" || species == "HO2"|| species == "H2O2" || species == "He"  || species == "Ne"
        energy = 0.0f0 * conversion_factor
    else
        throw(ArgumentError("Unknown species: $species"))
    end

    return energy
end

function get_electronic_energy_predis(species::String)

    """
    Return the electronic energy [J] for the given `species` for predissociation cases.
    """
    energy = 0.0f0
    conversion_factor = 1.602f-19

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
        throw(ArgumentError("Unknown species: $species"))
    end

    return energy
end

function allocate_velocity_new(reaction::PhotoReaction, E_excess, species_masses, p_photon)

    # Calculate velocities for the photodissociation products according to the moment and energy conservation quations
    # E_photon: Photon Energy in J
    # p_photon: linear momentum of photon in kg*m/s

    # 1. Get unitary vector for photon
    u_ph = reaction.sun_tuple

    # 2. Unpack species masses
    m_parent, m_heavy, m_light = species_masses

    # 3. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ m_parent .* reaction.v_parent

    # 4. Calculate unitary vectors for light and heavy products
    u_light = random_unit_tuple()
    u_heavy = random_unit_tuple()

    # 5. Calculate speed values for light and heavy products
    a = m_light * (m_heavy * (dot(u_heavy,u_ph)))^2 + (m_light^2 * m_heavy)*(dot(u_light, u_ph)^2)
    b = -2 * (dot(u_light, u_ph) * m_light * m_heavy) * (dot(total_momentum,u_ph))
    c = m_heavy * dot(total_momentum,u_ph)^2 - 2 * E_excess * (m_heavy * (dot(u_heavy,u_ph)))^2

    Δ = b^2 - 4 * a * c
    v_light_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_light_tuple = v_light_mag .* u_light

    v_heavy_tuple = (total_momentum .- m_light.*v_light_tuple) ./ m_heavy

    return v_light_tuple, v_heavy_tuple
end

# Simulate a single photodissociation reaction
function simulate_photodissociation(reaction::PhotoReaction, E_photon::Float32)

    # 1. Calculate photon linear momentum magnitude
    p_photon = calculate_photon_momentum(E_photon, reaction.sun_tuple)

    # 2. Get masses for all the species involved
    species_masses = get_masses(reaction.product_types[1], reaction.product_types[2], reaction.product_types[3])

    # 3. Calculate excess energy
    E_excess = calculate_excess_energy(reaction::PhotoReaction, species_masses[1], E_photon)

    # 4. Calculate velocity of heavy and light photoreaction product
    v_light_tuple, v_heavy_tuple = allocate_velocity_new(reaction, E_excess, species_masses, p_photon)

    return v_light_tuple, v_heavy_tuple
end

function simulate_photodissociation(reaction::PhotoReaction, E_photon::Real)
    return simulate_photodissociation(reaction, Float32(E_photon))
end

# Function to simulate multiple photodissociation reactions
function multiple_photodissociation(reaction::PhotoReaction, energy_vector::Vector{Float32})

    # This function is called when simulating a great number of times a specific photoreaction for a specific wavelength range to compare velocity values vs literature

    println("Simulating photodissociation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *""* " + " *""* reaction.product_names[3])

    final_speeds_light = []
    final_speeds_heavy = []

    # 1. Loop over photon energy vector
    for E_photon in energy_vector
        if E_photon > reaction.E_bond
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_light_tuple, v_heavy_tuple = simulate_photodissociation(reaction, E_photon)
            push!(final_speeds_light, v_light_tuple./1000)
            push!(final_speeds_heavy, v_heavy_tuple./1000)
        end
    end

    # 2. Show mean, STD and median speeds for photo products
    if reaction.display_info
        show_info(reaction, final_speeds_heavy, final_speeds_light)
    end

    return final_speeds_light, final_speeds_heavy
end

function multiple_photodissociation(reaction::PhotoReaction, energy_vector::AbstractVector{<:Real})
    return multiple_photodissociation(reaction, Float32.(collect(energy_vector)))
end

function show_info(reaction::PhotoReaction, final_speeds_heavy, final_speeds_light)
    final_speeds_heavy_norms = [norm(p) for p in final_speeds_heavy]
    final_speeds_light_norms = [norm(p) for p in final_speeds_light]

    data_speeds = DataFrame(
    Product = [reaction.product_names[2], reaction.product_names[3]],
    Mean_Speed = [mean(final_speeds_heavy_norms), mean(final_speeds_light_norms)],
    Median_Speed = [median(final_speeds_heavy_norms), median(final_speeds_light_norms)],
    STD_half = [std(final_speeds_heavy_norms)/2, std(final_speeds_light_norms)/2])
    println(data_speeds)
    println("")
end

end

