module SimplePhotoionisation

using Random, LinearAlgebra, Distributions, Statistics, DataFrames

const c = 299792458        # Speed of light in m/s
const h  = 6.62607015e-34  # Plank Constant in J.s
const m_el = 9.1093837e-31 # Electron mass in kg
const m_fund = 1.66054e-27 # 1 M.U.

struct PhotoReaction
    E_ionisation_eV::Float64 # Ionisation threshold energy in eV
    parent_velocity::Union{Tuple{Float64, Float64, Float64}, Float64} # Velocity of the parent molecule in m/s (H2O, OH, H2) -> USER INPUT
    sun_tuple::Union{Tuple{Float64, Float64, Float64}, Nothing}
    wvl_range::Union{Vector{Float64}, Nothing} # Contains photon wavelength for specified range. Not necessary if just running simple_photodissociation / simple_ionisation function
    energy_vector::Union{Vector{Float64}, Nothing} # Contains photon energies for specified range. Not necessary if just running simple_photodissociation / simple_ionisation function
    product_names::Union{Tuple{String, String, String}, Nothing} # USER INPUT -> Dict containing all involved species names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
    product_types::Union{Tuple{String, String, String}, Nothing}
    E_ionisation::Float64  # CALCULATED -> Ionisation energy in J
    display_info::Bool # Set true if you want to print photoproduct velocity analysis at the end

    function PhotoReaction(E_ionisation_eV, parent_velocity, sun_tuple, wvl_range, energy_vector, product_names, display_info)
        E_ionisation = E_ionisation_eV * 1.602e-19
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(E_ionisation_eV, parent_velocity, sun_tuple, wvl_range, energy_vector, product_names, product_types, E_ionisation, display_info)
    end
end

function random_unit_tuple()
    θ, φ = 2π * rand(), π * rand()
    return (sin(φ) * cos(θ), sin(φ) * sin(θ), cos(φ))
end

function get_masses(parent_name)
    # Get masses for involved photoreaction
    possible_parents = ("H", "H2", "OH", "H2O", "HO2", "H2O2", "He", "Ne")
    mass_dict = (1 * m_fund, 2 * m_fund, 17 * m_fund, 18 * m_fund, 33 * m_fund, 34 * m_fund, 4 * m_fund, 20 * m_fund)
    m_parent = m_ion = mass_dict[findfirst(isequal(parent_name), possible_parents)]

    return m_parent, m_ion
end

function allocate_velocity_new(reaction, E_photon, p_photon)

    # Calculate velocities for the photoionisation products according to the moment and energy conservation quations
    # E_photon: Photon Energy in J
    # p_photon: linear momentum of photon in kg*m/s
    # 1. Get masses for the parent species (also the mass of the product ion!)
    m_parent, m_ion = get_masses(reaction.product_types[1])

    # 2. Calculate excess energy for reaction
    E_excess = E_photon - reaction.E_ionisation
    
    # 3. Calculate unitary photon momentum tuple
    u_ph = p_photon ./ norm(p_photon)

    # 4. Generate random direction for parent molecule if not given and calculate velocity tuple
    if reaction.parent_velocity isa Float64
        v_parent = reaction.parent_velocity .* random_unit_tuple()
    elseif reaction.parent_velocity isa Tuple{Float64, Float64, Float64}
        v_parent = reaction.parent_velocity
    end

    # 5. Calculate total initial (conserved) momentum
    total_momentum = p_photon .+ m_parent .* v_parent

    # 6. Calculate unitary tuples for electron and product ion
    u_el = random_unit_tuple()
    u_ion = random_unit_tuple()

    # 7. Calculate speed values for electron and ion
    a = m_el * (m_ion * (dot(u_ion,u_ph)))^2 + (m_el^2 * m_ion)*(dot(u_el, u_ph)^2)
    b = -2 * (dot(u_el, u_ph) * m_el * m_ion) * (dot(total_momentum,u_ph))
    c = m_ion * dot(total_momentum,u_ph)^2 - 2 * E_excess * (m_ion * (dot(u_ion,u_ph)))^2

    Δ = b^2 - 4 * a * c
    v_el_mag = Δ ≥ 0 ? (-b + sqrt(Δ)) / (2 * a) : 0
    v_el_tuple= v_el_mag .* u_el

    v_ion_tuple= (total_momentum .- m_el.*v_el_tuple) ./ m_ion

    # 8. Only interested in ion speed

    return v_ion_tuple
end

# Simulate a single photodissociation reaction
function simulate_photoionisation(reaction::PhotoReaction, E_photon)

    # 1. Calculate photon linear momentum magnitude
    p_photon_magnitude = (E_photon) / c

    # 2. Calculate photon linear momentum vector by generating random incoming direction
    if reaction.sun_tuple isa Tuple{Float64, Float64, Float64}
        p_photon = p_photon_magnitude .* sun_tuple
    elseif reaction.sun_tuple isa Nothing
        p_photon = p_photon_magnitude .* random_unit_tuple()
    end

    # 3. Calculate velocity of product ion
    v_ion_tuple = allocate_velocity_new(reaction, E_photon, p_photon)

    # 5. Output trajectories and veloicties for product ion
    return v_ion_tuple
end

# Function to simulate multiple photodissociation reactions
function multiple_photoionisation(reaction::PhotoReaction)

    println("Simulating photoionisation reaction: " *""* reaction.product_names[1] *""* " + γ - > " *""* reaction.product_names[2] *" + e-")

    final_speeds_ion = []

    # 1. Loop over photon energy vector
    for (index, E_photon) in enumerate(reaction.energy_vector)

        if E_photon > reaction.E_ionisation
            # 1.1. Simulate individual photoreaction for every incoming photon
            v_ion_tuple= simulate_photoionisation(reaction, E_photon)

            # 1.2. Store trajectories and speeds of every ion product
            push!(final_speeds_ion, v_ion_tuple./1000)

        end
    end

    final_speeds_ion_norm = [norm(p) for p in final_speeds_ion]

    # 4. Show mean, STD and median speeds for product ions
    if reaction.display_info
        data_speeds = DataFrame(
        Product = [reaction.product_names[1]],
        Mean_Speed = [mean(final_speeds_ion_norm)],
        Median_Speed = [median(final_speeds_ion_norm)],
        STD_half = [std(final_speeds_ion_norm)/2])
        println(data_speeds)
        println("")
    end

    return final_speeds_ion
end

end