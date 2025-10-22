module solar_spectrum

using Interpolations

include(joinpath(@__DIR__, "solar_database.jl"))
using .solar_database
include("../../src/database/photodatabase.jl")
using .photodatabase

export flux_outputs

const c = 299792458 # Speed of light in m/s
const h  = 6.62607015f-34 # m2 kg / s


""" 
#:: FUNCTION: normalize_flux_distribution(wavelengths, solar_fluxes, wvl_range)
-------------------------------------------------------------------------------------------

INPUTS: 
- wavelengths::  solar_wavelength data from Huebner (2015)
- solar_fluxes::  Solar fluxes in photons/s (already includes cross section dependency)
- min_wavelength:: Lower bound for wavelength interval to study
- max_wavelength:: Upper bound for wavelength interval to study

OUTPUTS:
- filtered_wavelengths:: new wavelength tuple
- normalized_fluxes:: weighted probabilities for new wavelength vector
"""

function normalize_flux_distribution(wavelengths::NTuple{N, Int32}, solar_fluxes::NTuple{N, Float32}, wvl_range::NTuple{2, Float32})  where {N}

    @assert length(wavelengths) == length(solar_fluxes) "wavelengths and solar_fluxes must have the same length"

    # Convert tuples to vectors for indexing
    wl_vec = collect(wavelengths)
    flux_vec = collect(solar_fluxes)

    # Filter data to get wavelengths within the specified range
    filtered_indices = (wl_vec .>= wvl_range[1]) .& (wl_vec .<= wvl_range[2])
    filtered_wavelengths = wl_vec[filtered_indices]
    filtered_fluxes = flux_vec[filtered_indices]

    # Exclude Lyman-alpha if the condition holds
    if wvl_range[1] != 1.0 && wvl_range[2] != 95000.0
        if wvl_range[1] < 1208.0 && wvl_range[2] > 1220.0 # Filter out Lyman-alpha from 984-1357 A range!!
            filtered_wavelengths, filtered_fluxes = exclude_lyman_alpha(filtered_wavelengths, filtered_fluxes)
        end
    end

    # Normalize the fluxes
    normalized_fluxes = filtered_fluxes / sum(filtered_fluxes)

    return Tuple(filtered_wavelengths), Tuple(normalized_fluxes)
end

function normalize_flux_distribution(wavelengths::AbstractVector{<:Real}, solar_fluxes::AbstractVector{<:Real}, wvl_range::AbstractVector{<:Real})
    @assert length(wvl_range) == 2 "wvl_range should be of length 2!"
    N = length(wavelengths)
    return normalize_flux_distribution(ntuple(i -> Int32(wavelengths[i]), N), ntuple(i -> Float32(solar_fluxes[i]), N), ntuple(i -> Float32(wvl_range[i]), 2))
end

function normalize_flux_distribution(wavelengths::AbstractVector{<:Real}, solar_fluxes::AbstractVector{<:Real}, wvl_range::Tuple{<:Real,<:Real})
    N = length(wavelengths)
    return normalize_flux_distribution(ntuple(i -> Int32(wavelengths[i]), N), ntuple(i -> Float32(solar_fluxes[i]), N), (Float32(wvl_range[1]), Float32(wvl_range[2])))
end

function normalize_flux_distribution(wavelengths::NTuple{N, Real}, solar_fluxes::NTuple{N, Real}, wvl_range::Tuple{<:Real,<:Real}) where {N}
    return normalize_flux_distribution(Int32.(wavelengths), Float32.(solar_fluxes), (Float32(wvl_range[1]), Float32(wvl_range[2])))
end


function exclude_lyman_alpha(filtered_wavelengths::Vector{Int32}, filtered_fluxes::Vector{Float32})
    # Filter out Lyman-alpha from 984-1357 A range!!
    exclude_indices = findall(x -> x != 1212, filtered_wavelengths)  # Indices excluding 1212
    return filtered_wavelengths[exclude_indices], filtered_fluxes[exclude_indices]
end


"""
#:: FUNCTION: generate_weighted_random_wvl(wavelengths, normalized_fluxes)
-------------------------------------------------------------------------------------------

OBJECTIVE: Generate random wavelength according to the given weights derived from the species dependent photon flux in photons/s

INPUTS: 
- wavelengths::  wavelength tuple
- normalized_fluxes:: weighted probabilities for new wavelength vector

OUTPUTS:
- wavelength value in Armstrong
"""

function generate_weighted_random_wvl(wavelengths::NTuple{N,Int32}, normalized_fluxes::NTuple{N,Float32}) where {N}
    cumulative_distribution = Base.cumsum(collect(normalized_fluxes))
    random_value = rand(Float32)
    selected_index = findfirst(x -> x >= random_value, cumulative_distribution)
    return wavelengths[selected_index]
end

function generate_weighted_random_wvl(wavelengths::AbstractVector{<:Real}, normalized_fluxes::AbstractVector{<:Real})
    @assert length(wavelengths) == length(normalized_fluxes) "wavelengths and normalized_fluxes must have the same length"
    N = length(wavelengths)
    return generate_weighted_random_wvl(ntuple(i -> Int32(wavelengths[i]), N), ntuple(i -> Float32(normalized_fluxes[i]), N))
end


"""
#:: FUNCTION: generate_energy_vector(wavelengths, normalized_fluxes, num_reactions)
-------------------------------------------------------------------------------------------

OBJECTIVE: Generate vector of size num_reactions with random energies and wavelengths according to the given weights derived from the species dependent photon flux in photons/s

INPUTS: 
- wavelengths::  wavelength tuple
- normalized_fluxes:: weighted probabilities for new wavelength vector
- num_reactions:: number of photons to be simulated at once (usually 1)

OUTPUTS:
- wvl_vector, e_vector
"""

function generate_energy_vector(wavelengths::NTuple{N,Int32}, normalized_fluxes::NTuple{N,Float32}, num_reactions::Int) where {N}

    @assert 0 < num_reactions "Number of reactions must be at least 1!"

    if num_reactions == 1
        photon_wvl = generate_weighted_random_wvl(wavelengths, normalized_fluxes) # in Armstrong
        photon_energy = h * c / (photon_wvl*1f-10)
        return photon_wvl, photon_energy
    else
        wvl_vector = Float32[]
        e_vector = Float32[]
        for i in 1:num_reactions
            photon_wvl = generate_weighted_random_wvl(wavelengths, normalized_fluxes)  # in Armstrong
            push!(wvl_vector, photon_wvl)
            push!(e_vector, h * c / (photon_wvl*1f-10))
        end
        return wvl_vector, e_vector
    end
end


function generate_energy_vector(wavelengths::AbstractVector{<:Real}, normalized_fluxes::AbstractVector{<:Real}, num_reactions::Int)
    @assert length(wavelengths) == length(normalized_fluxes) "wavelengths and normalized_fluxes must have the same length"
    N = length(wavelengths)
    return generate_energy_vector(ntuple(i -> Int32(wavelengths[i]), N), ntuple(i -> Float32(normalized_fluxes[i]), N), num_reactions)
end


function get_solar_fluxes(parent_type::String, wvl_range::NTuple{2, Float64})

    @assert wvl_range[2] > wvl_range[1] "Wavelength bounds incorrect"

    if wvl_range == (1.0, 95000.0)
        wavelengths = solar_wavelength
        flux_quiet, flux_active = get_normalized_fluxes(parent_type)
    else
        flux_quiet_standard, flux_active_standard = get_standard_fluxes(parent_type)
        wavelengths, flux_quiet = normalize_flux_distribution(solar_wavelength, flux_quiet_standard, wvl_range)
        wavelengths, flux_active = normalize_flux_distribution(solar_wavelength, flux_active_standard, wvl_range)
    end

    return wavelengths, flux_quiet, flux_active
end

function get_solar_fluxes(parent_type::String, wvl_range::Tuple{<:Real,<:Real})
    return get_solar_fluxes(parent_type, Float32.(Tuple(wvl_range)))
end


""" 
#:: FUNCTION: flux_outputs(parent_type, wvl_range, solar_activity, num_reactions)
-------------------------------------------------------------------------------------------

OBJECTIVE: Get photon energy and wavelength for a certain parent species and degree of solar activity

INPUTS: 
- parent_type:: Parent species to be photodissociated/ionized
- wvl_range::  2D Tuple
- solar_activity:: decimal degree of solar activity (0.0 -> quiet sun, 1.0 -> active sun, or in between)
- num_reactions:: number of photons to be simulated at once (usually 1)

OUTPUTS:
- Final wavelength and energy vector / value
"""

function flux_outputs(parent_type::String, wvl_range::NTuple{2, Float64}, solar_activity::Float32, num_reactions::Int)

    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"

    wavelengths, flux_quiet, flux_active = get_solar_fluxes(parent_type, wvl_range)

    # Interpolate between quiet and active
    total_flux = (1.0f0 - solar_activity) .* flux_quiet .+ solar_activity .* flux_active

    # Return photon energy and wavelength vectors (size 1 if num_reactions is 1, vector otherwise)
    return generate_energy_vector(wavelengths, total_flux, num_reactions)

end

function flux_outputs(parent_type::String, wvl_range::AbstractVector{<:Real}, solar_activity::Real, num_reactions::Int)
    return flux_outputs(parent_type, Float64.(Tuple(wvl_range)), Float32(solar_activity), num_reactions)
end

function flux_outputs(parent_type::String, wvl_range::Tuple{<:Real,<:Real}, solar_activity::Real, num_reactions::Int)
    return flux_outputs(parent_type, Float64.(Tuple(wvl_range)), Float32(solar_activity), num_reactions)
end

end