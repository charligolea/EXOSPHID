# ─────────────────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────────────

""" 
    calculate_excess_energy_ionisation(E_bond, photon_energy)
-------------------------------------------------------------------------------------------
# Arguments
- `E_bond::Float32` -> Ionisation energy in J
- `photon_energy::Float32` -> In in J

# OBJECTIVE: 
- Calculate excess energy with a simplified approach
- For ionization, we neglect vibrorotational contributions compared to the energy of the 
    photon and electron

# Output:
- Excess energy in J
"""
function calculate_excess_energy_ionisation(E_bond::Float32, photon_energy::Float32)
    return photon_energy - E_bond
end



""" 
    allocate_velocity_ionisation(reaction, E_excess, species_masses, p_photon)
-------------------------------------------------------------------------------------------
# Arguments
- `reaction::PhotoReaction` -> PhotoReaction object
- `E_excess::Real` -> in J
- `species_masses::Float64` -> species scalar mass of parent species in kg
- `p_photon::NTuple{3, Real}` -> 3D Tuple in kg*m/s

# OBJECTIVE: 
- Calculate velocities for the photoionisation products according to the moment and energy 
    conservation equations
- The derivations of the quadratic equation can be consulted in the EXOSPHID WIKI

# Output:
- `v_ion_tuple`: 3D Tuple containing velocity components for the ionised parent
"""
function allocate_velocity_ionisation(reaction::PhotoReaction, E_excess::Real, 
                                    species_masses::Float64, p_photon::NTuple{3, Real})

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
    a = m_el * (m_ion * (dot(u_ion,u_ph)))^2 + 
        (m_el^2 * m_ion)*(dot(u_el, u_ph)^2)
    b = -2 * (dot(u_el, u_ph) * m_el * m_ion) * (dot(total_momentum,u_ph))
    c = m_ion * dot(total_momentum,u_ph)^2 - 
        2 * E_excess * (m_ion * (dot(u_ion,u_ph)))^2

    Δ = b^2 - 4 * a * c
    v_el_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_el_tuple= v_el_mag .* u_el

    v_ion_tuple= (total_momentum .- m_el.* v_el_tuple) ./ m_ion

    return Float32.(v_ion_tuple)
end


"""
    simulate_photoionisation(reaction, photon_energy)
-------------------------------------------------------------------------------------------
# Arguments
- `reaction::PhotoReaction`
- `photon_energy::Float32` -> In in J

# OBJECTIVE: 
- Simulate a single photoionisation reaction that has previously been determined from the 
    database

# Output: Output from `allocate_velocity()`
- `v_ion_tuple`: 3D Tuple containing velocity components for the ionised parent
"""
function simulate_photoionisation(reaction::PhotoReaction, photon_energy::Float32)

    # 1. Calculate photon linear momentum magnitude (kg*m/s)
    p_photon = calculate_photon_momentum(photon_energy, reaction.sun_tuple)

    # 2. Get masses for the parent species (also the mass of the product ion!)
    species_masses = get_masses(reaction.product_types[1]; mode="PI")

    # 3. Calculate excess energy for reaction (J)
    E_excess = calculate_excess_energy_ionisation(reaction.E_bond, photon_energy)

    # 4. Calculate velocity of product ion
    v_ion_tuple = 
        allocate_velocity_ionisation(reaction, E_excess, species_masses, p_photon)

    return v_ion_tuple
end

function simulate_photoionisation(reaction::PhotoReaction, photon_energy::Real)
    return simulate_photoionisation(reaction, Float32(photon_energy))
end


"""
    multiple_photoionisation(reaction, photon_energy_vector)
-------------------------------------------------------------------------------------------
# Arguments
- `reaction::PhotoReaction`
- `photon_energy_vector::Vector{Float32}` -> contains N scalar values of photon energy in J 
    in the wavelength range of choice

# OBJECTIVE: 
- Function to simulate multiple ionisation reactions
- Particularly interesting for validation studies ehere we want to generate multiple photons 
    at the same time for a specific parent 
AND wavelength range  AND reaction type and compare to literature values

# Output: Outputs from `allocate_velocity()`
- `final_speeds_ion`: Array of Size N. Every element is a 3D Tuple containing velocity 
    components for the ionised parent
"""
function multiple_photoionisation(reaction::PhotoReaction, 
                                photon_energy_vector::Vector{Float32})

    if reaction.display_info
        println("Simulating photoionisation reaction: " 
            *""* reaction.product_names[1] *""* " + γ - > " 
            *""* reaction.product_names[2] *" + e-")
    end

    final_speeds_ion = []

    # 1. Loop over photon energy vector
    for photon_energy in photon_energy_vector

        if photon_energy > reaction.E_bond
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_ion_tuple= simulate_photoionisation(reaction, photon_energy)

            # 1.2. Store trajectories and speeds of every ion product
            push!(final_speeds_ion, v_ion_tuple./1000)

        end
    end

    # 2. Show mean, STD and median speeds for product ions
    if reaction.display_info
        show_info_ionisation(reaction, final_speeds_ion)
    end

    return final_speeds_ion
end

function multiple_photoionisation(reaction::PhotoReaction, 
                                photon_energy_vector::AbstractVector{<:Real})
    return multiple_photoionisation(reaction, Float32.(collect(photon_energy_vector)))
end


"""
    show_info_ionisation(reaction, final_speeds_heavy, final_speeds_light)
-------------------------------------------------------------------------------------------
# Arguments
- `reaction::PhotoReaction`
- `final_speeds_ion::Vector{NTuple{3, Float64}}` -> output from `multiple_photoionisation`

# OBJECTIVE: 
- For the multiple ionisation case, show statistics of mean, median and STD speeds
- Only if `display_info` is set true
"""
function show_info_ionisation(reaction::PhotoReaction, final_speeds_ion::Vector{Any})
    final_speeds_ion_norm = norm.(final_speeds_ion)

    data_speeds = DataFrame(
    Product = [reaction.product_names[1]],
    Mean_Speed = [mean(final_speeds_ion_norm)],
    Median_Speed = [median(final_speeds_ion_norm)],
    STD_half = [std(final_speeds_ion_norm)/2])
    println(data_speeds)
    println("")
end


# ─────────────────────────────────────────────────────────────────────────────────────────
# EXPORTS
# ─────────────────────────────────────────────────────────────────────────────────────────

export simulate_photoionisation
export multiple_photoionisation