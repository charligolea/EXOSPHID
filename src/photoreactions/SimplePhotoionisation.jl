module SimplePhotoionisation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

const c = 2.99792458f8        # Speed of light in m/s
const m_el = 9.1093837e-31 # Electron mass in kg
const m_fund = 1.66054e-27 # 1 M.U.

struct PhotoReaction
    E_ionisation::Float32 # Ionisation threshold energy in J
    v_parent::NTuple{3, Float32} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::NTuple{3, Float32} # 
    product_names::NTuple{3, String} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::NTuple{3, String}
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_ionisation::Real, v_parent::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_ionisation), Float32.(v_parent), Float32.(sun_tuple), product_names, product_types, display_info)
    end

    function PhotoReaction(E_ionisation::Real, v_parent::Real, sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        vp = Float32(v_parent) .* random_unit_tuple()
        new(Float32(E_ionisation), vp, Float32.(sun_tuple), product_names, product_types, display_info)
    end
    
    function PhotoReaction(E_ionisation::Real, v_parent::NTuple{3, Float32}, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_ionisation), Float32.(v_parent), random_unit_tuple(), product_names, product_types, display_info)
    end

    function PhotoReaction(E_ionisation::Real, v_parent::Real, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        vp = Float32(v_parent) .* random_unit_tuple()
        new(Float32(E_ionisation), vp, random_unit_tuple(), product_names, product_types, display_info)
    end
end

function random_unit_tuple()
    θ, φ = Float32(2π) * rand(Float32), Float32(π) * rand(Float32)
    return (Float32(sin(φ) * cos(θ)), Float32(sin(φ) * sin(θ)), Float32(cos(φ)))
end

function get_masses(parent_name)
    # Get masses for involved photoreaction
    possible_parents = ("H", "H2", "OH", "H2O", "HO2", "H2O2", "He", "Ne")
    mass_dict = (1 * m_fund, 2 * m_fund, 17 * m_fund, 18 * m_fund, 33 * m_fund, 34 * m_fund, 4 * m_fund, 20 * m_fund)
    return mass_dict[findfirst(isequal(parent_name), possible_parents)]
end

calculate_photon_momentum(E_photon::Float32, sun_tuple) = (E_photon/c) .* sun_tuple
calculate_excess_energy(reaction::PhotoReaction, E_photon::Float32) = E_photon - reaction.E_ionisation

function allocate_velocity_new(reaction, E_excess, species_masses, p_photon)

    # Calculate velocities for the photoionisation products according to the moment and energy conservation quations
    # E_photon: Photon Energy in J
    # p_photon: linear momentum of photon in kg*m/s

    # 1. Get unitary vector for photon
    u_ph = reaction.sun_tuple

    # 2. Unpack species masses
    m_parent = m_ion = species_masses

    # 3. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ m_parent .* reaction.v_parent

    # 4. Calculate unitary tuples for electron and product ion
    u_el = random_unit_tuple()
    u_ion = random_unit_tuple()

    # 5. Calculate speed values for electron and ion
    a = m_el * (m_ion * (dot(u_ion,u_ph)))^2 + (m_el^2 * m_ion)*(dot(u_el, u_ph)^2)
    b = -2 * (dot(u_el, u_ph) * m_el * m_ion) * (dot(total_momentum,u_ph))
    c = m_ion * dot(total_momentum,u_ph)^2 - 2 * E_excess * (m_ion * (dot(u_ion,u_ph)))^2

    Δ = b^2 - 4 * a * c
    v_el_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_el_tuple= v_el_mag .* u_el

    v_ion_tuple= (total_momentum .- m_el.*v_el_tuple) ./ m_ion

    return v_ion_tuple
end

# Simulate a single photodissociation reaction
function simulate_photoionisation(reaction::PhotoReaction, E_photon::Float32)

    # 1. Calculate photon linear momentum magnitude
    p_photon = calculate_photon_momentum(E_photon, reaction.sun_tuple)

    # 2. Get masses for the parent species (also the mass of the product ion!)
    species_masses = get_masses(reaction.product_types[1])

    # 3. Calculate excess energy for reaction
    E_excess = calculate_excess_energy(reaction, E_photon)

    # 2. Calculate velocity of product ion
    v_ion_tuple = allocate_velocity_new(reaction, E_excess, species_masses, p_photon)

    return v_ion_tuple
end

function simulate_photoionisation(reaction::PhotoReaction, E_photon::Real)
    return simulate_photoionisation(reaction, Float32(E_photon))
end

# Function to simulate multiple photodissociation reactions
function multiple_photoionisation(reaction::PhotoReaction, energy_vector::Vector{Float32})

    println("Simulating photoionisation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *" + e-")

    final_speeds_ion = []

    # 1. Loop over photon energy vector
    for E_photon in energy_vector

        if E_photon > reaction.E_ionisation
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_ion_tuple= simulate_photoionisation(reaction, E_photon)

            # 1.2. Store trajectories and speeds of every ion product
            push!(final_speeds_ion, v_ion_tuple./1000)

        end
    end

    # 2. Show mean, STD and median speeds for product ions
    if reaction.display_info
        show_info(reaction, final_speeds_ion)
    end

    return final_speeds_ion
end

function multiple_photoionisation(reaction::PhotoReaction, energy_vector::AbstractVector{<:Real})
    return multiple_photoionisation(reaction, Float32.(collect(energy_vector)))
end

function show_info(reaction::PhotoReaction, final_speeds_ion)
    final_speeds_ion_norm = [norm(p) for p in final_speeds_ion]

    data_speeds = DataFrame(
    Product = [reaction.product_names[1]],
    Mean_Speed = [mean(final_speeds_ion_norm)],
    Median_Speed = [median(final_speeds_ion_norm)],
    STD_half = [std(final_speeds_ion_norm)/2])
    println(data_speeds)
    println("")
end

end