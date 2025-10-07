module SimplePhotodissociation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

const c = 299792458        # Speed of light in m/s
const h  = 6.62607015e-34  # Plank Constant in J.s
const m_fund = 1.66054e-27 # 1 M.U.

struct PhotoReaction
    E_bond_eV::Float64 # Threshold energy for given photodissociation reaction in eV -> USER INPUT
    parent_velocity::Union{Tuple{Float64, Float64, Float64}, Float64} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::Union{Tuple{Float64, Float64, Float64}, Nothing} # 
    wvl_range::Union{Vector{Float64}, Nothing} # Contains photon wavelength for specified range. Not necessary if just running simple_photodissociation / simple_ionisation function
    energy_vector::Union{Vector{Float64}, Nothing} # Contains photon energies for specified range. Not necessary if just running simple_photodissociation / simple_ionisation function
    product_names::Union{Tuple{String, String, String}, Nothing} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::Union{Tuple{String, String, String}, Nothing}
    E_bond::Float64 # CALCULATED -> Photon Energy in J
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_bond_eV, parent_velocity, sun_tuple, wvl_range, energy_vector, product_names, display_info)
        E_bond = E_bond_eV * 1.602e-19
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(E_bond_eV, parent_velocity, sun_tuple, wvl_range, energy_vector, product_names, product_types, E_bond, display_info)
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
    elseif species == "OH(X2π)" || species == "OH"
        energy = 0.20 * conversion_factor
    elseif species == "OH(A2Σ+)"
        energy = 0.25 * conversion_factor
    elseif species == "OH(1Σ+)" || species == "OH(12Δ/22Π)"
        energy = 0.25 * conversion_factor
    elseif species == "OH(A2Σ+, v'=3)"
        # energy = 1.50 * conversion_factor
        energy = 0.00 * conversion_factor # TO REVIEW
    elseif species == "OH(A2Σ+, v'=2)"
        energy = 1.05 * conversion_factor
        # energy = 1.10 * conversion_factor
    elseif species == "OH(B2Σ)"
        energy = 0.05 * conversion_factor
    elseif species == "OH(D2Σ)"
        energy = 0.20 * conversion_factor
    elseif species == "O" || species == "O(3P)" || species == "O(1D)" || species == "O(1S)"
        energy = 0.0 * conversion_factor
    elseif species == "H2"
        energy = 0.30 * conversion_factor
    elseif species == "H(1s)" || species == "H(2s,2p)" || species == "H"  || species == "He"  || species == "Ne"
        energy = 0.0 * conversion_factor
    elseif species == "H(-)"
        energy = 0.0 * conversion_factor
    elseif species == "HO2"
        energy = 0.0 * conversion_factor
    elseif species == "H2O2"
        energy = 0.0 * conversion_factor
    else
        throw(ArgumentError("Unknown species: $species"))
    end

    return energy
end

function get_electronic_energy_predis(species::String)

    """
    Return the electronic energy [J] for the given `species` for predissociation cases.
    """
    energy = 0.0
    conversion_factor = 1.602e-19

    if species == "OH(X2π)" || species == "OH"
        energy = 0.00 * conversion_factor
    elseif species == "OH(A2Σ+)" || species == "OH(A2Σ+, v'=3)" || species == "OH(A2Σ+, v'=2)"
        energy = 4.05 * conversion_factor
    elseif species == "OH(1Σ+)"
        energy = 4.05 * conversion_factor
    elseif species == "OH(12Δ/22Π)"
        energy = 6.50 * conversion_factor
    elseif species == "OH(B2Σ)"
        energy = 8.65 * conversion_factor
    elseif species == "OH(D2Σ)"
        energy = 10.18 * conversion_factor
    else
        throw(ArgumentError("Unknown species: $species"))
    end

    return energy
end


function allocate_velocity_new(reaction, E_photon, p_photon)

    # Calculate velocities for the photodissociation products according to the moment and energy conservation quations
    # E_photon: Photon Energy in J
    # p_photon: linear momentum of photon in kg*m/s

    # 1. Get masses for all the species involved
    m_parent, m_heavy, m_light = get_masses(reaction.product_types[1], reaction.product_types[2], reaction.product_types[3])

    # 2. Calculate excess energy for reaction    
    if reaction.product_types[1] == "OH"
        E_parent = 1/2 * m_parent * norm(reaction.parent_velocity)^2
        E_excess = (E_parent + E_photon - get_electronic_energy_predis(reaction.product_names[1]) - get_vibrorotational_energy(reaction.product_names[1]))
    else
        E_parent = 1/2 * m_parent * norm(reaction.parent_velocity)^2 + get_vibrorotational_energy(reaction.product_names[1])
        E_excess = (E_photon + E_parent) - reaction.E_bond - (get_vibrorotational_energy(reaction.product_names[2]) + get_vibrorotational_energy(reaction.product_names[3]))
    end

    # 3. Calculate unitary photon momentum tuple
    u_ph = p_photon ./ norm(p_photon)

    # 4. Calculate parent velocity tuple
    if reaction.parent_velocity isa Float64
        v_parent = reaction.parent_velocity .* random_unit_tuple()
    elseif reaction.parent_velocity isa Tuple{Float64, Float64, Float64}
        v_parent = reaction.parent_velocity
    end

    # 5. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ m_parent .* v_parent

    # 6. Calculate unitary vectors for light and heavy products
    u_light = random_unit_tuple()
    u_heavy = random_unit_tuple()

    # 7. Calculate speed values for light and heavy products
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
function simulate_photodissociation(reaction::PhotoReaction, E_photon)

    # 1. Calculate photon linear momentum magnitude
    p_photon_magnitude = (E_photon) / c

    # 2. Calculate photon linear momentum vector by generating random incoming direction
    if reaction.sun_tuple isa Tuple{Float64, Float64, Float64}
        p_photon = p_photon_magnitude .* sun_tuple
    elseif reaction.sun_tuple isa Nothing
        p_photon = p_photon_magnitude .* random_unit_tuple()
    end

    # 3. Calculate velocity of heavy and light photoreaction product
    v_light_tuple, v_heavy_tuple = allocate_velocity_new(reaction, E_photon, p_photon)

    # 5. Output trajectories and veloicties for "light" and "heavy" photo product
    return v_light_tuple, v_heavy_tuple
end

# Function to simulate multiple photodissociation reactions
function multiple_photodissociation(reaction::PhotoReaction)

    # This function is called when simulating a great number of times a specific photoreaction for a specific wavelength range to compare velocity values vs literature

    println("Simulating photodissociation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *""* " + " *""* reaction.product_names[3])

    final_speeds_light = []
    final_speeds_heavy = []

    # 1. Loop over photon energy vector

    for (index, E_photon) in enumerate(reaction.energy_vector)

        if E_photon > reaction.E_bond
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_light_tuple, v_heavy_tuple = simulate_photodissociation(reaction, E_photon)
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

