# ─────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────────────

""" 
    calculate_excess_energy_dissociation(reaction, m_parent, photon_energy)
-------------------------------------------------------------------------------------------

# Arguments 
- `reaction::PhotoReaction`
- `m_parent:: Float64` -> mass of parent species in kg
- `photon_energy::Real` -> in J

# OBJECTIVE: 
- Calculate excess energy taking into consideration vibrorotational energies of parents and 
    products.
- OH is a special case, as it is predissociated by the photon. To understand the theory 
    behind predissociation and simulation chocies for EXOSPHID, refer to WIKI

# Output: 
- Excess energy in J
"""
function calculate_excess_energy_dissociation(reaction::PhotoReaction, m_parent::Float64, 
                                            photon_energy::Real)

    # PREDISSOCIATION
    if reaction.product_types[1] == "OH" && reaction.product_names[1] != "OH(DPD)"
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2
        E_excess = 
            E_parent + 
            photon_energy - 
            get_electronic_energy_predis(reaction.product_names[1]) - 
            get_vibrorotational_energy(reaction.product_names[1])

    # DOUBLE DISSOCIATION
    elseif occursin("(DPD)", reaction.product_names[1])
        E_parent = 0.5f0 * m_parent * norm(reaction.v_parent)^2
        E_excess = (photon_energy + E_parent) - reaction.E_bond

    # STANDARD CASE
    else
        E_parent = 
            0.5f0 * m_parent * norm(reaction.v_parent)^2 + 
            get_vibrorotational_energy(reaction.product_names[1])
        E_excess = 
            photon_energy + 
            E_parent - 
            reaction.E_bond -
            get_vibrorotational_energy(reaction.product_names[2]) - 
            get_vibrorotational_energy(reaction.product_names[3])

    end

    return E_excess
end


""" 
    allocate_velocity(reaction, E_excess, species_masses, p_photon)
-------------------------------------------------------------------------------------------

# Arguments
- `reaction::PhotoReaction`
- `E_excess::Real` -> in J
- `species_masses::NTuple{3, Float64}` -> 3D Tuple containing: 
    * [1]: parent mass
    * [2]: heavy product mass
    * [3]: light product mass
- `p_photon::NTuple{3, Real}` -> 3D Tuple in kg*m/s

# OBJECTIVE: 
- Calculate velocities for the photodissociation products according to the moment and 
    energy conservation equations
- The derivations of the quadratic equation can be consulted in the EXOSPHID WIKI

# Output:
- `v_heavy_tuple`: 3D Tuple containing velocity components for the heavier photolysis product 
    (e.g., for H2O -> OH + H, it would be OH)
- `v_light_tuple`: 3D Tuple containing velocity components for the lighter photolysis product 
    (e.g., for H2O -> OH + H, it would be H)
"""
function allocate_velocity_dissociation(reaction::PhotoReaction, 
                                    E_excess::Real, 
                                    species_masses::NTuple{3, Float64}, 
                                    p_photon::NTuple{3, Real})

    # 1. Get unitary vector for photon
    u_ph = reaction.sun_tuple

    # 2. Unpack species masses
    m_parent, m_heavy, m_light = species_masses

    # 3. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ m_parent .* reaction.v_parent

    # 4. Calculate unitary vectors for light and heavy products
    u_light = random_unit_tuple()
    u_light = dot(u_light, u_ph) ≥ 0 ? u_light : (-u_light[1], -u_light[2], -u_light[3])

    u_heavy = random_unit_tuple()
    u_heavy = dot(u_heavy, u_ph) ≥ 0 ? u_heavy : (-u_heavy[1], -u_heavy[2], -u_heavy[3])

    # 5. Calculate speed values for light and heavy products
    a = m_light * (m_heavy * (dot(u_heavy,u_ph)))^2 + 
        (m_light^2 * m_heavy)*(dot(u_light, u_ph)^2)
    b = -2 * (dot(u_light, u_ph) * m_light * m_heavy) * (dot(total_momentum,u_ph))
    c = m_heavy * dot(total_momentum,u_ph)^2 - 
        2 * E_excess * (m_heavy * (dot(u_heavy,u_ph)))^2

    Δ = b^2 - 4 * a * c
    v_light_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_light_tuple = v_light_mag .* u_light

    v_heavy_tuple = (total_momentum .- m_light.* v_light_tuple) ./ m_heavy

    return Float32.(v_light_tuple), Float32.(v_heavy_tuple)
end


"""
    simulate_photodissociation(reaction, photon_energy)
-------------------------------------------------------------------------------------------

# Arguments
- `reaction::PhotoReaction` object
- `photon_energy::Float32` -> in J

# OBJECTIVE: 
- Simulate a single photodissociation reaction that has previously been determined from the 
    database

# Output: Outputs from `allocate_velocity()`
- `v_heavy_tuple`: 3D Tuple containing velocity components for the heavier photolysis product 
    (e.g., for H2O -> OH + H, it would be OH)
- `v_light_tuple`: 3D Tuple containing velocity components for the lighter photolysis product 
    (e.g., for H2O -> OH + H, it would be H)
"""
function simulate_photodissociation(reaction::PhotoReaction, photon_energy::Float32)

    # 1. Calculate photon linear momentum magnitude (kg*m/s)
    p_photon = calculate_photon_momentum(photon_energy, reaction.sun_tuple)

    # 2. Get masses for all the species involved
    species_masses = (get_masses(reaction.product_types[1]), 
                      get_masses(reaction.product_types[2]), 
                      get_masses(reaction.product_types[3]))

    # 3. Calculate excess energy (J)
    E_excess = calculate_excess_energy_dissociation(reaction::PhotoReaction, 
                                                species_masses[1], photon_energy)

    # 4. Calculate velocity of heavy and light photoreaction product
    v_light_tuple, v_heavy_tuple = 
        allocate_velocity_dissociation(reaction, E_excess, species_masses, p_photon)

    return v_light_tuple, v_heavy_tuple
end

function simulate_photodissociation(reaction::PhotoReaction, photon_energy::Real)
    return simulate_photodissociation(reaction, Float32(photon_energy))
end


"""
    multiple_photodissociation(reaction, photon_energy_vector)
-------------------------------------------------------------------------------------------

# Arguments
- `reaction::PhotoReaction` object
- `photon_energy_vector::Vector{Float32}` -> contains N scalar values of photon energy in J 
    in the wavelength range of choice

# OBJECTIVE: 
- Function to simulate multiple photodissociation reactions
- Particularly interesting for validation studies ehere we want to generate multiple photons 
    at the same time for a specific parent 
AND wavelength range  AND reaction type and compare to literature values

# Output: Outputs from `allocate_velocity()`
- `final_speeds_light`: Array of Size N. Every element is a 3D Tuple containing velocity 
    components for the heavier photolysis product (e.g., for H2O -> OH + H, it would be OH)
- `final_speeds_heavy`: Array of Size N. Every element is a 3D Tuple containing velocity 
    components for the lighter photolysis product (e.g., for H2O -> OH + H, it would be H)
"""
function multiple_photodissociation(reaction::PhotoReaction, 
                                photon_energy_vector::Vector{Float32})

    final_speeds_light = NTuple{3, Float32}[]
    final_speeds_heavy = NTuple{3, Float32}[]

    # 1. Loop over photon energy vector
    for photon_energy in photon_energy_vector
        if photon_energy > reaction.E_bond
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_light_tuple, v_heavy_tuple = 
                simulate_photodissociation(reaction, photon_energy)
            push!(final_speeds_light, v_light_tuple./1000)
            push!(final_speeds_heavy, v_heavy_tuple./1000)
        end
    end

    # 2. Show mean, STD and median speeds for photo products
    if reaction.display_info
        show_info_dissociation(reaction, final_speeds_heavy, final_speeds_light)
    end

    return final_speeds_light, final_speeds_heavy
end

function multiple_photodissociation(reaction::PhotoReaction, 
                                photon_energy_vector::AbstractVector{<:Real})
    return multiple_photodissociation(reaction, Float32.(collect(photon_energy_vector)))
end


"""
    show_info(reaction, final_speeds_heavy, final_speeds_light)
-------------------------------------------------------------------------------------------

# Arguments
- `reaction::PhotoReaction` object
- `final_speeds_heavy::Vector{Any}, final_speeds_light::Vector{Any}` -> outputs from 
    `multiple_photodissociation`

# OBJECTIVE: 
- For the multiple photodssociation case, show statistics of mean, median and STD speeds
- Only if display_info is set true
"""
function show_info_dissociation(reaction::PhotoReaction, 
                            final_speeds_heavy::Vector{Tuple{Float32, Float32, Float32}}, 
                            final_speeds_light::Vector{Tuple{Float32, Float32, Float32}})

    println("Simulating photodissociation reaction: " 
            *""* reaction.product_names[1] *""* " + γ - > " 
            *""* reaction.product_names[2] *""* " + " *""* reaction.product_names[3])

    final_speeds_heavy_norms = norm.(final_speeds_heavy)
    final_speeds_light_norms = norm.(final_speeds_light)

    data_speeds = DataFrame(
    Product = [reaction.product_names[2], reaction.product_names[3]],
    Mean_Speed = [mean(final_speeds_heavy_norms), mean(final_speeds_light_norms)],
    Median_Speed = [median(final_speeds_heavy_norms), median(final_speeds_light_norms)],
    STD_half = [std(final_speeds_heavy_norms)/2, std(final_speeds_light_norms)/2])
    println(data_speeds)
    println("")
end


# ─────────────────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ─────────────────────────────────────────────────────────────────────────────────────────

export simulate_photodissociation
export multiple_photodissociation