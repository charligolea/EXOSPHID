if !isdefined(Main, :EXOSPHID)
    include(joinpath(@__DIR__, "..", "src", "EXOSPHID.jl"))
    using .EXOSPHID
end

using StatsPlots
using Random
using LinearAlgebra
using KernelDensity
using StatsBase
using Plots


"""
    organize_product_velocities(prod_types, prod_vels)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Organize `photodestruction` runs into a Dict object associating every possible 
    photoproduct with its velocity tuples

# INPUTS:
- `prod_types`::Vector{Any} -> stores all the product_type arrays from running N times
    the `photodestruction` function
- `prod_vels``::Vector{Any} -> stores all the product_velocity arrays from running N times
    the `photodestruction` function
"""
function organize_product_velocities(prod_types::Vector{Any}, prod_vels::Vector{Any})

    speed_results = Dict{String, Vector{NTuple{3, Float32}}}()

    for (products, velocities) in zip(prod_types, prod_vels)
        @assert length(products) == length(velocities)

        for (species, vel) in zip(products, velocities)
            if !haskey(speed_results, species)
                speed_results[species] = NTuple{3,Float32}[]
            end
            push!(speed_results[species], vel./1000)
        end
    end
    
    return speed_results
end


"""
    plot_speed_distributions(speed_results, parent_type)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Plot velocities of photodestruction products and store in appropriate parent folder.

# INPUTS:
- `speed_results::Dict{String, Vector{NTuple{3, Float32}}}` -> Every key must correspond to
    a certain photoproduct, and its value is a vector containing all the associated velocity
    tuples
- `parent_type::String` -> Name of parent species. Must be in `exosphid_species``
"""
function plot_speed_distributions(speed_results::Dict{String, Vector{NTuple{3, Float32}}}, parent_type::String)
        
    Random.seed!(42)

    color_dict = Dict(k => palette(:viridis, length(speed_results))[i]
                    for (i, k) in enumerate(keys(speed_results)))

    for (species, vels) in speed_results
        isempty(vels) && continue
        speeds = norm.(vels)

        # Remove 1–99% outliers
        lb, ub = quantile(speeds, [0.01, 0.99])
        filtered = filter(s -> lb ≤ s ≤ ub, speeds)
        length(filtered) < 3 && continue

        # Histogram
        h = fit(Histogram, filtered, nbins=30)
        w = h.weights ./ sum(h.weights)
        centers = @. 0.5*(h.edges[1][1:end-1] + h.edges[1][2:end])
        widths = diff(h.edges[1])

        p = bar(centers, w, width=widths, label="", legend=false,
                xlabel="$species Speed (km/s)", ylabel="Relative Frequency",
                fillcolor=color_dict[species], fillalpha=0.4, linecolor=:transparent)

        # KDE scaled to histogram height
        xs = range(minimum(filtered), maximum(filtered), length=500)
        ys = pdf(kde(filtered), xs)
        plot!(p, xs, ys .* (maximum(w) / maximum(ys)), linewidth=2,
            linecolor=color_dict[species], label="")

        output_dir = joinpath(@__DIR__, "..", "..", "output", "speed_distributions", parent_type)
        mkpath(output_dir)

        savefig(p, joinpath(output_dir, "$species.pdf"))    
    end
end


"""
    speed_distributions_example_usage(parent_type, num_reactions)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Plot speed distributions of all the photodestruction products of a given parent species.
- Assumes high dt to make sure photon-parent interaction always occurs

# INPUTS:
- `parent_type::String` -> Name of parent species. Must be in `exosphid_species``
- `num_reactions::Integer` -> Number of photoon-parent interactions to simulate
"""
function speed_distributions_example_usage(parent_type::String, num_reactions::Integer)

    prod_types, prod_vels = [], []

    for i in 1:num_reactions
        ro, rn, pt, pv, wr = photodestruction(0, 1e10, parent_type, EXOSPHID.velocities[parent_type], nothing)
        if ro; push!(prod_types, pt); push!(prod_vels, pv); end
    end

    speed_results = organize_product_velocities(prod_types, prod_vels)
    plot_speed_distributions(speed_results, parent_type)
end

export speed_distributions_example_usage