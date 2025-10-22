module SimplePhotoionisation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

const c = 2.99792458f8        # Speed of light in m/s
const m_el = 9.1093837e-31 # Electron mass in kg

include("../database/photodatabase.jl")
using .photodatabase

struct PhotoReaction
    E_ionisation::Float32 # Ionisation threshold energy in J
    v_parent::NTuple{3, Float32} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::NTuple{3, Float32} # Solar vector === direction of incoming photons
    product_names::NTuple{3, String} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::NTuple{3, String} # Extracts (+) from strings for ionic species to get the "clean" string for the atomic species (e.g. H(+) -> H)
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_ionisation::Real, v_parent::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_ionisation), Float32.(v_parent), Float32.(sun_tuple), product_names, product_types, display_info)
    end

    function PhotoReaction(E_ionisation::Real, v_parent::Real, sun_tuple::NTuple{3, Real}, product_names::NTuple{3, String}, display_info::Bool)
        vp = Float32(v_parent) .* random_unit_tuple()
        PhotoReaction(Float32(E_ionisation), v_parent, Float32.(sun_tuple), product_names, display_info)
    end
    
    function PhotoReaction(E_ionisation::Real, v_parent::NTuple{3, Real}, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        PhotoReaction(Float32(E_ionisation), Float32.(v_parent), random_unit_tuple, product_names, display_info)
    end

    function PhotoReaction(E_ionisation::Real, v_parent::Real, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        vp = Float32(v_parent) .* random_unit_tuple()
        PhotoReaction(Float32(E_ionisation), vp, random_unit_tuple(), product_names, display_info)
    end

        function PhotoReaction(E_ionisation::Real, v_parent::AbstractArray{<:Real}, sun_tuple::AbstractArray{<:Real}, product_names::NTuple{3, String}, display_info::Bool)
        @assert length(v_parent) == 3 "v_parent must have length 3"
        @assert length(sun_tuple) == 3 "sun_tuple must have length 3"
        vp = Float32.(v_parent)
        st = Float32.(sun_tuple)
        PhotoReaction(Float32(E_ionisation), vp, st, product_names, display_info)
    end

    function PhotoReaction(E_ionisation::Real, v_parent::Real, sun_tuple::AbstractArray{<:Real}, product_names::NTuple{3, String}, display_info::Bool)
        vp = Float32(v_parent) .* random_unit_tuple()
        @assert length(sun_tuple) == 3 "sun_tuple must have length 3"
        st = Float32.(sun_tuple)
        PhotoReaction(Float32(E_ionisation), vp, st, product_names, display_info)
    end

    function PhotoReaction(E_ionisation::Real, v_parent::AbstractArray{<:Real}, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        @assert length(v_parent) == 3 "v_parent must have length 3"
        vp = Float32.(v_parent)
        st = random_unit_tuple()
        PhotoReaction(Float32(E_ionisation), vp, st, product_names, display_info)
    end

end

"""
#:: FUNCTION: random_unit_tuple()

# OBJECTIVE: For the cases where parent velocity has been provided as scalar, or solar vector has not been provided, generate random unitary vector
"""

function random_unit_tuple()
    θ, φ = Float32(2π) * rand(Float32), Float32(π) * rand(Float32)
    return (Float32(sin(φ) * cos(θ)), Float32(sin(φ) * sin(θ)), Float32(cos(φ)))
end


""" 
#:: FUNCTION: calculate_photon_momentum(E_photon, sun_tuple)
#-------------------------------------------------------------------------------------------
# Arguments
- E_photon: in J  

# Output: Momentum vector for the photon in kg*m/s
"""

calculate_photon_momentum(E_photon::Float32, sun_tuple) = (E_photon/c) .* sun_tuple

""" 
#:: FUNCTION: calculate_excess_energy(E_ionisation, E_photon)
#-------------------------------------------------------------------------------------------
# Arguments
- E_ionisation:: Ionisation energy in J
- E_photon in J

# OBJECTIVE: 
- Calculate excess energy with a simplified approach
- For ionization, we neglect vibrorotational contributions compared to the energy of the photon and electron

# Output: Excess energy in J
"""

calculate_excess_energy(E_ionisation::Float32, E_photon::Float32) = E_photon - E_ionisation



""" 
#:: FUNCTION: allocate_velocity(reaction, E_excess, species_masses, p_photon)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- E_excess in J
- species_masses: scalar mass of parent species in kg
- p_photon: 3D Tuple in kg*m/s

# OBJECTIVE: 
- Calculate velocities for the photoionisation products according to the moment and energy conservation equations
- The derivations of the quadratic equation can be consulted in the EXOSPHID WIKI

# Output:
- v_ion_tuple: 3D Tuple containing velocity components for the ionised parent
"""

function allocate_velocity(reaction::PhotoReaction, E_excess, species_masses::Float64, p_photon)

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

    return Float32.(v_ion_tuple)
end


"""
#:: FUNCTION: simulate_photoionisation(reaction, E_photon)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- E_photon in J

# OBJECTIVE: 
- Simulate a single photoionisation reaction that has previously been determined from the database

# Output: Output from allocate_velocity()
- v_ion_tuple: 3D Tuple containing velocity components for the ionised parent
"""

function simulate_photoionisation(reaction::PhotoReaction, E_photon::Float32)

    # 1. Calculate photon linear momentum magnitude (kg*m/s)
    p_photon = calculate_photon_momentum(E_photon, reaction.sun_tuple)

    # 2. Get masses for the parent species (also the mass of the product ion!)
    species_masses = get_masses(reaction.product_types[1]; mode="PI")

    # 3. Calculate excess energy for reaction (J)
    E_excess = calculate_excess_energy(reaction.E_ionisation, E_photon)

    # 4. Calculate velocity of product ion
    v_ion_tuple = allocate_velocity(reaction, E_excess, species_masses, p_photon)

    return v_ion_tuple
end

function simulate_photoionisation(reaction::PhotoReaction, E_photon::Real)
    return simulate_photoionisation(reaction, Float32(E_photon))
end


"""
#:: FUNCTION: multiple_photoionisation(reaction, energy_vector)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- energy_vector: contains N scalar values of photon energy in J in the wavelength range of choice

# OBJECTIVE: 
- Function to simulate multiple ionisation reactions
- Particularly interesting for validation studies ehere we want to generate multiple photons at the same time for a specific parent 
AND wavelength range  AND reaction type and compare to literature values

# Output: Outputs from allocate_velocity()
- final_speeds_ion: Array of Size N. Every element is a 3D Tuple containing velocity components for the ionised parent
"""

function multiple_photoionisation(reaction::PhotoReaction, energy_vector::Vector{Float32})

    if reaction.display_info
        println("Simulating photoionisation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *" + e-")
    end

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


"""
#:: FUNCTION: show_info(reaction, final_speeds_heavy, final_speeds_light)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- final_speeds_ion: output from multiple_photoionisation

# OBJECTIVE: 
- For the multiple ionisation case, show statistics of mean, median and STD speeds
- Only if display_info is set true
"""

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