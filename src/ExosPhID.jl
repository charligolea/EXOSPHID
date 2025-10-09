module ExosPhID

using Random, LinearAlgebra, Distributions

include("photoreactions/SimplePhotodissociation.jl")
using .SimplePhotodissociation

include("photoreactions/SimplePhotoionisation.jl")
using .SimplePhotoionisation

include("database/solar_spectrum.jl")
using .solar_spectrum

include("database/photodatabase.jl")
include("photoreactions/photodestruction_logic.jl")

export photodestruction

""" 
# Arguments
- solar_activity = From 0 (Quiet Sun) to 1 (Active Sun)
- dt: analyzed simulation time window
- parent_type:: String with parent molecule type ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2")
- parent_velocity:: Tuple{Float32, Float32, Float32} or Float32, Parent velocity vector in m/s    

# Outputs
- reaction_occurence::Boolean
- present_reaction:: String poiting at the reaction that is occuring
- product_types:: ["T1", "T2", "T3"]
- product_velocities:: [(vx1, vy1, vz1), (vx2, vy2, vz2), (vx3, vy3, vz3)]
"""

function random_unit_tuple()
    θ, φ = 2f0 * Float32(π) * rand(Float32), Float32(π) * rand(Float32)
    return (Float32(sin(φ) * cos(θ)), Float32(sin(φ) * sin(θ)), Float32(cos(φ)))
end

function photodestruction(solar_activity::Float32, dt::Float32, parent_type::String, parent_velocity::Tuple{Float32, Float32, Float32}, sun_tuple::Tuple{Float32, Float32, Float32})

    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"
    @assert 0.0 <= dt "dt must be positive!"
    @assert parent_type in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne") "Invalid parent species: $parent_type"

    # 0. Check if reaction is happening according to probability
    k = get_photodestruction_rates(parent_type, solar_activity) # Photodestruction rate
    reaction_occurence = is_photoreaction_occuring(k, dt)

    if reaction_occurence == true

        # 1. Generate Photon Energy (J) and Wavelength (A) from Huebner Solar Flux Distributions (1992, 2015)
        photon_wvl, photon_energy = flux_outputs(parent_type, (1.0f0, 95000.0f0), solar_activity, Int32(1))

        # 2. Simulate photoreaction
        wvl_threshold = get_wvl_threshold(parent_type)

        if photon_wvl <= wvl_threshold
            photochemical_info = get_species_photochemical_info(parent_type)
            photoreaction_characteristics = get_photoreaction_characteristics(photon_wvl, photochemical_info)

            if  photoreaction_characteristics[1] !== nothing
                product_velocities, product_types = call_photodestruction_logic(photoreaction_characteristics, photochemical_info, parent_velocity, sun_tuple, photon_energy)
                reaction_name = photoreaction_characteristics[2]
                wvl_range = photoreaction_characteristics[4]
                return reaction_occurence, reaction_name, product_types, product_velocities, wvl_range
            else 
                return false, "", [], [], ()
            end
        else
            return false, "", [], [], ()
        end
    else
        return false, "", [], [], ()
    end

end


function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Real, sun_tuple::Tuple)
    @assert 0.0 <= parent_velocity "Parent velocity must be positive!"
    pv = Float32(parent_velocity) .* random_unit_tuple()
    st = Tuple(Float32.(sun_tuple))
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Real, sun_tuple::Nothing)
    @assert 0.0 <= parent_velocity "Parent velocity must be positive!"
    pv = Float32(parent_velocity) .* random_unit_tuple()
    st = random_unit_tuple()
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Tuple, sun_tuple::Nothing)
    pv = Tuple(Float32.(parent_velocity))
    st = random_unit_tuple()
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Tuple, sun_tuple::Tuple)
    pv = Tuple(Float32.(parent_velocity))
    st = Tuple(Float32.(sun_tuple))
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

end