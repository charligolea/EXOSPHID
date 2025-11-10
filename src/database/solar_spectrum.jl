""" 
`h`: Planck Constant in m^2 kg / s
"""
const h  = 6.62607015f-34 # m2 kg / s

""" 
    normalize_flux_distribution(wavelengths, solar_fluxes, wvl_range)
-------------------------------------------------------------------------------------------

# INPUTS: 
- `wavelengths::NTuple{N, Integer}` ->  `solar_wavelength` data from Huebner (2015)
- `solar_fluxes::NTuple{N, Float32}` ->  Solar fluxes in photons/s (already includes cross 
    section dependency)
- `wvl_range::NTuple{2, Float32})` -> Lower and upper bound for wavelength range, expressed 
    in Angstrom

# OUTPUTS:
- `filtered_wavelengths` -> new wavelength tuple
- `normalized_fluxes` -> weighted probabilities for new wavelength vector
"""
function normalize_flux_distribution(wavelengths::NTuple{N, Integer}, 
                                    solar_fluxes::NTuple{N, Float32}, 
                                    wvl_range::NTuple{2, Float32})  where {N}

    @assert length(wavelengths) == length(solar_fluxes) "wavelengths and solar_fluxes must
                                                        have the same length"

    # Convert tuples to vectors for indexing
    wl_vec = collect(wavelengths)
    flux_vec = collect(solar_fluxes)

    # Filter data to get wavelengths within the specified range
    filtered_indices = (wl_vec .>= wvl_range[1]) .& (wl_vec .<= wvl_range[2])
    filtered_wavelengths = wl_vec[filtered_indices]
    filtered_fluxes = flux_vec[filtered_indices]

    # Exclude Lyman-alpha if the condition holds
    if wvl_range[1] != 1.0 && wvl_range[2] != 95000.0
        # Filter out Lyman-alpha from 984-1357 A range!!
        if wvl_range[1] < 1208.0 && wvl_range[2] > 1220.0 
            filtered_wavelengths, filtered_fluxes = 
                exclude_lyman_alpha(filtered_wavelengths, filtered_fluxes)
        end
    end

    # Normalize the fluxes
    normalized_fluxes = filtered_fluxes / sum(filtered_fluxes)

    return Tuple(filtered_wavelengths), Tuple(normalized_fluxes)
end

function normalize_flux_distribution(wavelengths::AbstractVector{<:Real}, 
                                    solar_fluxes::AbstractVector{<:Real}, 
                                    wvl_range::AbstractVector{<:Real})
    @assert length(wvl_range) == 2 "wvl_range should be of length 2!"
    N = length(wavelengths)
    return normalize_flux_distribution(ntuple(i -> Int(wavelengths[i]), N), 
        ntuple(i -> Float32(solar_fluxes[i]), N), 
        ntuple(i -> Float32(wvl_range[i]), 2))
end

function normalize_flux_distribution(wavelengths::AbstractVector{<:Real}, 
                                    solar_fluxes::AbstractVector{<:Real}, 
                                    wvl_range::Tuple{<:Real,<:Real})
    N = length(wavelengths)
    return normalize_flux_distribution(ntuple(i -> Int(wavelengths[i]), N), 
        ntuple(i -> Float32(solar_fluxes[i]), N), 
        (Float32(wvl_range[1]), Float32(wvl_range[2])))
end

function normalize_flux_distribution(wavelengths::NTuple{N, Real},
                                    solar_fluxes::NTuple{N, Real}, 
                                    wvl_range::Tuple{<:Real,<:Real}) where {N}
    return normalize_flux_distribution(Int.(wavelengths), 
        Float32.(solar_fluxes), 
        (Float32(wvl_range[1]), Float32(wvl_range[2])))
end


function exclude_lyman_alpha(filtered_wavelengths::Vector{Integer}, 
                            filtered_fluxes::Vector{Float32})
    # Filter out Lyman-alpha from 984-1357 A range!!
    exclude_indices = findall(x -> x != 1212, filtered_wavelengths)  
    return filtered_wavelengths[exclude_indices], filtered_fluxes[exclude_indices]
end


"""
    generate_weighted_random_wvl(wavelengths, normalized_fluxes)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Generate random wavelength according to the given weights derived from the species 
    dependent photon flux in photons/s

# INPUTS: 
- `wavelengths::NTuple{N,Integer}` ->  wavelength tuple
- `normalized_fluxes::NTuple{N,Float32})` -> weighted probabilities for new wavelength 
    vector

OUTPUTS:
- wavelength value in Armstrong
"""
function generate_weighted_random_wvl(wavelengths::NTuple{N,Integer}, 
                                    normalized_fluxes::NTuple{N,Float32}) where {N}
    cumulative_distribution = Base.cumsum(collect(normalized_fluxes))
    random_value = rand(Float32)
    selected_index = findfirst(x -> x >= random_value, cumulative_distribution)
    return wavelengths[selected_index]
end

function generate_weighted_random_wvl(wavelengths::AbstractVector{<:Real}, 
                                    normalized_fluxes::AbstractVector{<:Real})
    @assert length(wavelengths) == length(normalized_fluxes) "wavelengths and 
                                normalized_fluxes must have the same length"
    L = length(wavelengths)
    return generate_weighted_random_wvl(ntuple(i -> Int(wavelengths[i]), L), 
        ntuple(i -> Float32(normalized_fluxes[i]), L))
end

function generate_weighted_random_wvl(wavelengths::AbstractVector{<:Real}, 
                                    normalized_fluxes::NTuple{N, Real}) where {N}
    @assert length(wavelengths) == length(normalized_fluxes) "wavelengths and 
                                normalized_fluxes must have the same length"
    L = length(wavelengths)
    return generate_weighted_random_wvl(ntuple(i -> Int(wavelengths[i]), L), 
        Float32.(normalized_fluxes))
end

function generate_weighted_random_wvl(wavelengths::NTuple{N, Real}, 
                                    normalized_fluxes::AbstractVector{<:Real}) where {N}
    @assert length(wavelengths) == length(normalized_fluxes) "wavelengths and 
                                normalized_fluxes must have the same length"
    L = length(wavelengths)
    return generate_weighted_random_wvl(Int.(wavelengths), 
        ntuple(i -> Float32(normalized_fluxes[i]), L))
end


"""
    generate_energy_vector(wavelengths, normalized_fluxes, num_reactions)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- # Generate vector of size num_reactions with random energies and wavelengths according to 
    the given weights derived from the species dependent photon flux in photons/s

# INPUTS: 
- `wavelengths::NTuple{N,Integer}` ->  wavelength tuple
- `normalized_fluxes::NTuple{N,Float32}` -> weighted probabilities for new wavelength vector
- `num_reactions::Integer` -> number of photons to be simulated at once (usually 1)

# OUTPUTS:
If the number of reactions is 1
- `wvl_vector` -> vector containing wavelength values in Angstrom
- `photon_energy_vector` -> vector containing photon energies in J

If the number of reactions is larger than 1:
- `photon_wvl` -> Photon Wavelength in Angstrom
- `photon_energy` -> Photon Energy in J

"""
function generate_energy_vector(wavelengths::NTuple{N,Integer}, 
                                normalized_fluxes::NTuple{N,Float32}, 
                                num_reactions::Integer) where {N}

    @assert 0 < num_reactions "Number of reactions must be at least 1!"

    if num_reactions == 1
        photon_wvl = generate_weighted_random_wvl(wavelengths, normalized_fluxes) # in A
        photon_energy = h * c / (photon_wvl*1f-10)
        return photon_wvl, photon_energy
    else
        wvl_vector = Float32[]
        photon_energy_vector = Float32[]
        for i in 1:num_reactions
            photon_wvl = generate_weighted_random_wvl(wavelengths, normalized_fluxes) # in A
            push!(wvl_vector, photon_wvl)
            push!(photon_energy_vector, h * c / (photon_wvl*1f-10))
        end
        return wvl_vector, photon_energy_vector
    end
end


function generate_energy_vector(wavelengths::AbstractVector{<:Real}, 
                                normalized_fluxes::AbstractVector{<:Real}, 
                                num_reactions::Integer)
    @assert length(wavelengths) == length(normalized_fluxes) "wavelengths and 
                                normalized_fluxes must have the same length"
    N = length(wavelengths)
    return generate_energy_vector(ntuple(i -> Int(wavelengths[i]), N), 
        ntuple(i -> Float32(normalized_fluxes[i]), N), 
        num_reactions)
end

function generate_energy_vector(wavelengths::AbstractVector{<:Real}, 
                                normalized_fluxes::Tuple{Vararg{<:Real}}, 
                                num_reactions::Integer)
    N = length(normalized_fluxes)
    @assert length(wavelengths) == N "wavelengths and normalized_fluxes 
                                    must have the same length"
    return generate_energy_vector(ntuple(i -> Int(wavelengths[i]), N), 
        ntuple(i -> Float32(normalized_fluxes[i]), N), 
        num_reactions)
end


"""
    get_solar_fluxes(parent_type, wvl_range)
-------------------------------------------------------------------------------------------

# OBJECTIVE: Get normalized fluxes for a given parent species and ina specific wavelength 
    range

# INPUTS: 
- `wvl_range::NTuple{2, Real}` ->  2D Tuple
- `parent_type::String` -> within `exosphid_species`

# OUTPUTS:
- `wvl_vector` -> vector containing wavelength values in Angstrom
- `e_vector` -> vector containing photon energies in J
"""

function get_solar_fluxes(parent_type::String, wvl_range::NTuple{2, Real})

    @assert wvl_range[2] > wvl_range[1] "Wavelength bounds incorrect"

    if wvl_range == (1.0, 95000.0)
        wavelengths = solar_wavelength
        flux_quiet, flux_active = get_normalized_fluxes(parent_type)
    else
        assert_not_same_interval(wvl_range, solar_wavelength)
        flux_quiet_standard, flux_active_standard = get_standard_fluxes(parent_type)
        wavelengths, flux_quiet = 
            normalize_flux_distribution(solar_wavelength, flux_quiet_standard, wvl_range)
        wavelengths, flux_active = 
            normalize_flux_distribution(solar_wavelength, flux_active_standard, wvl_range)
    end

    return wavelengths, flux_quiet, flux_active
end

function get_solar_fluxes(parent_type::String, wvl_range::AbstractVector{<:Real})
    @assert len(wvl_range) == 2 "Wavelength range array must be of size 2"
    return get_solar_fluxes(parent_type, Tuple(wvl_range))
end

"""
    assert_not_same_interval(wrange, solar_wavelength)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Checks that upper and lower wavlength bounds are not exactly between 2 consecutive 
        wavelength bins, otherwise the generated arrays would be of size 0.

# INPUTS:
- `wrange`:: Tuple{<:Real,<:Real} -> wavelength range bounds
- `solar_wavelength`::Tuple -> reference wavelength array
"""
function assert_not_same_interval(wrange::Tuple{<:Real,<:Real}, solar_wavelength::Tuple)
    solar_vec = collect(solar_wavelength)
    idxs = map(w -> searchsortedfirst(solar_vec, w), wrange)

    @assert !(all(i -> 1 < i <= length(solar_vec), idxs) && 
        idxs[1] == idxs[2] && 
        (solar_vec[idxs[1]] != wrange[1] && 
        solar_vec[idxs[2]] != wrange[2]))
        "Upper and lower wavelength range $wrange lie within 
        2 consecutive solar_wavelength bins"
end


""" 
    flux_outputs(parent_type, wvl_range, solar_activity, num_reactions)
-------------------------------------------------------------------------------------------

# OBJECTIVE: Get photon energy and wavelength for a certain parent species and degree of 
    solar activity

# INPUTS: 
- `parent_type::String` -> Parent species to be photodissociated/ionized
- `wvl_range::NTuple{2, Real}` ->  2D Tuple
- `solar_activity::Float32` -> decimal degree of solar activity 
    (0.0 -> quiet sun, 1.0 -> active sun, or in between)
- `num_reactions::Int` -> number of photons to be simulated at once (usually 1)

# OUTPUTS:
- Final wavelength and energy vector / value
"""
function flux_outputs(parent_type::String, wvl_range::NTuple{2, Real}, 
                    solar_activity::Float32, num_reactions::Int)

    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"

    wavelengths, flux_quiet, flux_active = get_solar_fluxes(parent_type, wvl_range)

    # Interpolate between quiet and active
    total_flux = (1.0f0 - solar_activity) .* flux_quiet .+ solar_activity .* flux_active

    # Return photon energy and wavelength vectors 
    # (size 1 if num_reactions is 1, vector otherwise)
    return generate_energy_vector(wavelengths, total_flux, num_reactions)

end

function flux_outputs(parent_type::String, wvl_range::AbstractVector{<:Real}, 
                    solar_activity::Real, num_reactions::Int)
    return flux_outputs(parent_type, Tuple(wvl_range), 
        Float32(solar_activity), num_reactions)
end

function flux_outputs(parent_type::String, wvl_range::Tuple{<:Real,<:Real}, 
                    solar_activity::Real, num_reactions::Int)
    return flux_outputs(parent_type, Tuple(wvl_range), 
        Float32(solar_activity), num_reactions)
end


export flux_outputs
export get_solar_fluxes