""" 
#:: FUNCTION: photodestruction(solar_activity, dt, parent_type, parent_velocity, sun_tuple)
#-------------------------------------------------------------------------------------------
# Arguments
- solar_activity = From 0 (Quiet Sun) to 1 (Active Sun)
- dt: analyzed simulation time window
- parent_type:: String with parent molecule type. See exosphid_species variable in photodatabase.jl for possible species
- parent_velocity:: 3D array or Scalar, Parent velocity vector in m/s    

# Outputs
- reaction_occurence::Boolean
- reaction_name:: Identifier for the reaction (e.g. H2O-PD1 is the first PhotoDissociation mechanism for H2O )
- product_types:: ["T1", "T2", "T3"]
- product_velocities:: [(vx1, vy1, vz1), (vx2, vy2, vz2), (vx3, vy3, vz3)]
"""

"""
#:: FUNCTION: random_unit_tuple()

# OBJECTIVE: For the cases where parent velocity has been provided as scalar, or solar vector has not been provided, generate random unitary vector
"""
function random_unit_tuple()
    θ, φ = 2f0 * Float32(π) * rand(Float32), Float32(π) * rand(Float32)
    return (Float32(sin(φ) * cos(θ)), Float32(sin(φ) * sin(θ)), Float32(cos(φ)))
end

function photodestruction(solar_activity::Float32, dt::Float32, parent_type::String, parent_velocity::Tuple{Float32, Float32, Float32}, sun_tuple::Tuple{Float32, Float32, Float32})

    # 1. Get photochemical information for parent species
    photochemical_info = get_species_photochemical_info(parent_type)

    # 2. Check if reaction is happening according to probability
    k = get_photodestruction_rates(photochemical_info, solar_activity) # Photodestruction rate
    reaction_occurence = is_photoreaction_occuring(k, dt)

    if reaction_occurence == true

        # 3. Generate Photon Energy (J) and Wavelength (A) from Huebner Solar Flux Distributions (1992, 2015)
        photon_wvl, photon_energy = flux_outputs(parent_type, (1, 95000), solar_activity, 1)

        # 4. Check that photon is sufficiently energetic to trigger photolysis
        wvl_threshold = get_wvl_threshold(photochemical_info)

        if photon_wvl <= wvl_threshold
            # 5. Determine current reaction for given photon
            current_reaction = get_current_reaction(photon_wvl, photochemical_info)

            if  current_reaction.present_reaction !== nothing
                # 6. Simulate photodestruction process
                product_velocities, product_types = call_photodestruction_logic(current_reaction, photochemical_info, parent_velocity, sun_tuple, photon_energy)
                return reaction_occurence, current_reaction.reaction_name, product_types, product_velocities, current_reaction.wvl_range
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


# ─────────────────────────────────────────────
# 2. SCALAR VELOCITY + TUPLE SUN
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Real, sun_tuple::Tuple)
    @assert 0.0 <= parent_velocity "Parent velocity must be positive!"
    @assert len(sun_tuple) == 3 "Give solar vector as a 3D object! (preferably tuple)"
    pv = Float32(parent_velocity) .* random_unit_tuple()
    st = Tuple(Float32.(sun_tuple))
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

# ─────────────────────────────────────────────
# 3. SCALAR VELOCITY + NO SUN VECTOR
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Real, sun_tuple::Nothing)
    @assert 0.0 <= parent_velocity "Parent velocity must be positive!"
    pv = Float32(parent_velocity) .* random_unit_tuple()
    st = random_unit_tuple()
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

# ─────────────────────────────────────────────
# 4. TUPLE VELOCITY + NO SUN VECTOR
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Tuple, sun_tuple::Nothing)
    @assert len(parent_velocity) == 3 "Give parent_velocity as a 3D object or a scalar"
    pv = Tuple(Float32.(parent_velocity))
    st = random_unit_tuple()
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

# ─────────────────────────────────────────────
# 5. TUPLE VELOCITY + TUPLE SUN VECTOR
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Tuple, sun_tuple::Tuple)
    @assert len(parent_velocity) == 3 "Give parent_velocity as a 3D object or a scalar"
    @assert len(sun_tuple) == 3 "Give solar vector as a 3D object! (preferably tuple)"
    pv = Tuple(Float32.(parent_velocity))
    st = Tuple(Float32.(sun_tuple))
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

# ─────────────────────────────────────────────
# 6. REAL VELOCITY, VECTOR SUN TUPLE
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::Real, sun_tuple::AbstractVector{<:Real})
    @assert 0.0 <= parent_velocity "Parent velocity must be positive!"
    @assert len(sun_tuple) == 3 "Solar vector must be 3D object!"
    pv = Float32(parent_velocity) .* random_unit_tuple()
    st = Tuple(Float32.(sun_tuple))
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end


# ─────────────────────────────────────────────
# 7. VECTOR VELOCITY + VECTOR SUN
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::AbstractVector{<:Real}, sun_tuple::AbstractVector{<:Real})
    @assert length(parent_velocity) == 3 "parent_velocity must have length 3"
    @assert length(sun_tuple) == 3 "sun_tuple must have length 3"
    pv = Tuple(Float32.(parent_velocity))
    st = Tuple(Float32.(sun_tuple))
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end


# ─────────────────────────────────────────────
# 8. VECTOR VELOCITY + NO SUN VECTOR
# ─────────────────────────────────────────────
function photodestruction(solar_activity::Real, dt::Real, parent_type::String, parent_velocity::AbstractVector{<:Real}, sun_tuple::Nothing)
    @assert length(parent_velocity) == 3 "parent_velocity must have length 3"
    pv = Tuple(Float32.(parent_velocity))
    st = random_unit_tuple()
    return photodestruction(Float32(solar_activity), Float32(dt), parent_type, pv, st)
end

export photodestruction
export random_unit_tuple