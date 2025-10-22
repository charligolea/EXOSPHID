module SimplePhotodissociation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

include("../database/photodatabase.jl")
using .photodatabase

const c = 2.99792458f8        # Speed of light in m/s

struct PhotoReaction
    E_bond::Float32 # Threshold energy for given photodissociation reaction in J -> USER INPUT
    v_parent::NTuple{3, Float32} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::NTuple{3, Float32} # Solar vector === direction of incoming photons
    product_names::NTuple{3, String} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::NTuple{3, String} # Extracts elements in parenthesis, duch as electronic states, from product names (e.g. OH(X^2Pi) -> OH)
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_bond::Real, v_parent::NTuple{3, Float32}, sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_bond), Float32.(v_parent), Float32.(sun_tuple), product_names, product_types, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::Real, sun_tuple::NTuple{3, Real}, product_names::NTuple{3, String}, display_info::Bool)
        vp = Float32(v_parent) .* random_unit_tuple()
        PhotoReaction(Float32(E_bond), vp, Float32.(sun_tuple), product_names, display_info)
    end
    
    function PhotoReaction(E_bond::Real, v_parent::NTuple{3, Real}, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        st = random_unit_tuple()
        PhotoReaction(Float32(E_bond), Float32.(v_parent), st, product_names, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::Real, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        vp = Float32(v_parent) .* random_unit_tuple()
        st = random_unit_tuple()
        PhotoReaction(Float32(E_bond), vp, st, product_names, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::AbstractArray{<:Real}, sun_tuple::AbstractArray{<:Real}, product_names::NTuple{3, String}, display_info::Bool)
        @assert length(v_parent) == 3 "v_parent must have length 3"
        @assert length(sun_tuple) == 3 "sun_tuple must have length 3"
        vp = Float32.(v_parent)
        st = Float32.(sun_tuple)
        PhotoReaction(Float32(E_bond), vp, st, product_names, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::Real, sun_tuple::AbstractArray{<:Real}, product_names::NTuple{3, String}, display_info::Bool)
        vp = Float32(v_parent) .* random_unit_tuple()
        @assert length(sun_tuple) == 3 "sun_tuple must have length 3"
        st = Float32.(sun_tuple)
        PhotoReaction(Float32(E_bond), vp, st, product_names, display_info)
    end

    function PhotoReaction(E_bond::Real, v_parent::AbstractArray{<:Real}, sun_tuple::Nothing, product_names::NTuple{3, String}, display_info::Bool)
        @assert length(v_parent) == 3 "v_parent must have length 3"
        vp = Float32.(v_parent)
        st = random_unit_tuple()
        PhotoReaction(Float32(E_bond), vp, st, product_names, display_info)
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

calculate_photon_momentum(E_photon::Real, sun_tuple::NTuple{3, Float32}) = (E_photon/c) .* sun_tuple


""" 
#:: FUNCTION: calculate_excess_energy(reaction, m_parent, E_photon)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- m_parent: mass of parent species in kg
- E_photon in J

# OBJECTIVE: 
- Calculate excess energy taking into consideration vibrorotational energies of parents and products.
- OH is a special case, as it is predissociated by the photon. To understand the theory behind 
predissociation and simulation chocies for EXOSPHID, refer to WIKI

# Output: Excess energy in J
"""

function calculate_excess_energy(reaction::PhotoReaction, m_parent::Float64, E_photon::Real)

    if reaction.product_types[1] == "OH" && reaction.product_names[1] != "OH(DPD)"
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2
        E_excess = (E_parent + E_photon - get_electronic_energy_predis(reaction.product_names[1]) - get_vibrorotational_energy(reaction.product_names[1]))
    elseif occursin("(DPD)", reaction.product_names[1])
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2 + get_vibrorotational_energy(reaction.product_names[1])
        E_excess = (E_photon + E_parent) - reaction.E_bond 
    else
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2 + get_vibrorotational_energy(reaction.product_names[1])
        E_excess = (E_photon + E_parent) - reaction.E_bond - (get_vibrorotational_energy(reaction.product_names[2]) + get_vibrorotational_energy(reaction.product_names[3]))
    end

    return E_excess
end


""" 
#:: FUNCTION: allocate_velocity(reaction, E_excess, species_masses, p_photon)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- E_excess in J
- species_masses: 3D Tuple containing: [1]: parent mass, [2]: heavy product mass, [3]: light product mass
- p_photon: 3D Tuple in kg*m/s

# OBJECTIVE: 
- Calculate velocities for the photodissociation products according to the moment and energy conservation equations
- The derivations of the quadratic equation can be consulted in the EXOSPHID WIKI

# Output:
- v_heavy_tuple: 3D Tuple containing velocity components for the heavier photolysis product (e.g., for H2O -> OH + H, it would be OH)
- v_light_tuple: 3D Tuple containing velocity components for the lighter photolysis product (e.g., for H2O -> OH + H, it would be H)
"""

function allocate_velocity(reaction::PhotoReaction, E_excess, species_masses, p_photon)

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

    return Float32.(v_light_tuple), Float32.(v_heavy_tuple)
end


"""
#:: FUNCTION: simulate_photodissociation(reaction, E_photon)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- E_photon in J

# OBJECTIVE: 
- Simulate a single photodissociation reaction that has previously been determined from the database

# Output: Outputs from allocate_velocity()
- v_heavy_tuple: 3D Tuple containing velocity components for the heavier photolysis product (e.g., for H2O -> OH + H, it would be OH)
- v_light_tuple: 3D Tuple containing velocity components for the lighter photolysis product (e.g., for H2O -> OH + H, it would be H)
"""

function simulate_photodissociation(reaction::PhotoReaction, E_photon::Float32)

    # 1. Calculate photon linear momentum magnitude (kg*m/s)
    p_photon = calculate_photon_momentum(E_photon, reaction.sun_tuple)

    # 2. Get masses for all the species involved
    species_masses = get_masses(reaction.product_types[1]; heavy_child_name= reaction.product_types[2], light_child_name = reaction.product_types[3], mode="PD")

    # 3. Calculate excess energy (J)
    E_excess = calculate_excess_energy(reaction::PhotoReaction, species_masses[1], E_photon)

    # 4. Calculate velocity of heavy and light photoreaction product
    v_light_tuple, v_heavy_tuple = allocate_velocity(reaction, E_excess, species_masses, p_photon)

    return v_light_tuple, v_heavy_tuple
end

function simulate_photodissociation(reaction::PhotoReaction, E_photon::Real)
    return simulate_photodissociation(reaction, Float32(E_photon))
end


"""
#:: FUNCTION: simulate_photodissociation(reaction, energy_vector)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- energy_vector: contains N scalar values of photon energy in J in the wavelength range of choice

# OBJECTIVE: 
- Function to simulate multiple photodissociation reactions
- Particularly interesting for validation studies ehere we want to generate multiple photons at the same time for a specific parent 
AND wavelength range  AND reaction type and compare to literature values

# Output: Outputs from allocate_velocity()
- final_speeds_light: Array of Size N. Every element is a 3D Tuple containing velocity components for the heavier photolysis product (e.g., for H2O -> OH + H, it would be OH)
- final_speeds_heavy: Array of Size N. Every element is a 3D Tuple containing velocity components for the lighter photolysis product (e.g., for H2O -> OH + H, it would be H)
"""

function multiple_photodissociation(reaction::PhotoReaction, energy_vector::Vector{Float32})

    if reaction.display_info
        println("Simulating photodissociation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *""* " + " *""* reaction.product_names[3])
    end

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


"""
#:: FUNCTION: show_info(reaction, final_speeds_heavy, final_speeds_light)
#-------------------------------------------------------------------------------------------
# Arguments
- reaction:: PhotoReaction object
- final_speeds_heavy, final_speeds_light: outputs from multiple_photodissociation

# OBJECTIVE: 
- For the multiple photodssociation case, show statistics of mean, median and STD speeds
- Only if display_info is set true
"""

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

