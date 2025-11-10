# ─────────────────────────────────────────────────────────────────────────────────────────
# Main Struct for Ionisation and Dissociation Reactions
# ─────────────────────────────────────────────────────────────────────────────────────────

"""
    PhotoReaction()

- `E_bond::Float32` -> Threshold energy for given photodissociation or photoionisation 
    reaction in J -> USER INPUT
- `v_parent::NTuple{3, Float32}` -> Velocity of the parent molecule in m/s (H2O, OH, H2)
- `sun_tuple::NTuple{3, Float32}` -> Solar vector === direction of incoming photons
- `product_names::NTuple{3, String}` -> USER INPUT -> Dict containing all involved species 
    names. Must contain 3 keys: "parent_name", "heavy_child_name", "light_child_name""
- `product_types::NTuple{3, String}` -> Extracts elements in parenthesis, duch as electronic 
    states, from product names (e.g. OH(X^2Pi) -> OH)
- `display_info::Bool` -> Set true if you want to print photoproduct velocity analysis
"""
struct PhotoReaction
    E_bond::Float32
    v_parent::NTuple{3, Float32}
    sun_tuple::NTuple{3, Float32}
    product_names::NTuple{3, String}
    product_types::NTuple{3, String}
    display_info::Bool

    function PhotoReaction(E_bond::Real, v_parent::NTuple{3, Float32}, 
            sun_tuple::NTuple{3, Float32}, product_names::NTuple{3, String}, 
            display_info::Bool)
        product_types = map(s -> replace(s, r"\(.*\)" => ""),product_names)
        new(Float32(E_bond), Float32.(v_parent), 
            Float32.(sun_tuple), product_names, product_types, 
            display_info)
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
    calculate_photon_momentum(E_photon, sun_tuple)
-------------------------------------------------------------------------------------------
# Arguments
- `E_photon::Real` -> in J  
- `sun_tuple::NTuple{3, Float32}` -> Solar vector === direction of incoming photons

# Output:
- Momentum vector for the photon in kg*m/s
"""
function calculate_photon_momentum(E_photon::Real, sun_tuple::NTuple{3, Float32})
    return (E_photon/c) .* sun_tuple
end