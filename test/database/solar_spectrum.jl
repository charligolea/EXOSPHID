using Random

println("TESTING solar_spectrum.jl ...............")


""" 
    sample_two_separated(vals)
-------------------------------------------------------------------------------------------

# OBJECTIVE:
- Generate a 2D tuple with minimum and maximum wavelength values to avoid landing exactly
   between 2 consecutive solar_wavelength elements

# OUTPUT:
- `wvl_range::NTuple{2, Float32})` -> Lower and upper bound for wavelength range, expressed 
    in Angstrom
"""
function sample_two_separated()
    vals = sort(collect(EXOSPHID.solar_wavelength))

    # define intervals between consecutive values
    intervals = [(vals[i], vals[i+1]) for i in 1:length(vals)-1]

    # pick two different intervals
    chosen_intervals = randperm(length(intervals))[1:2]

    # sample one number in each interval
    n1 = rand() * (intervals[chosen_intervals[1]][2] - intervals[chosen_intervals[1]][1]) + 
            intervals[chosen_intervals[1]][1]
    n2 = rand() * (intervals[chosen_intervals[2]][2] - intervals[chosen_intervals[2]][1]) + 
            intervals[chosen_intervals[2]][1]

    return (Float32(min(n1, n2)), Float32(max(n1, n2)))
end

@testset verbose=true "solar_spectrum.jl" begin

    # Even if you select a smaller range of wvl check that sorted and same size

    for pt in exosphid_species

        @testset "$pt" verbose=true begin 
                
            @testset "Check flux Array Lengths are correct for reduced WVL range" begin
                for i in 1:1000 # Try out 1000 cases
                    wvrange = sample_two_separated()
                    tol = 1e-5
                    wvl, fq, fa = get_solar_fluxes(pt, wvrange)

                    @test length(fq) == length(wvl)
                    @test length(fa) == length(wvl)
                end
            end

            @testset "Check Normalization is well done for reduced WVL range" begin
                for i in 1:1000 # Try out 1000 cases
                    wvrange = sample_two_separated()
                    tol = 1e-5
                    wvl, fq, fa = get_solar_fluxes(pt, wvrange)

                    @test abs(sum(fq) - 1.0) < tol
                    @test abs(sum(fa) - 1.0) < tol
                end
            end 
        end
    end
end


println("............... COMPLETED TESTING solar_spectrum.jl\n")