println("TESTING of solar_database.jl ...............  ")

@testset verbose=true "solar_database.jl" begin


    @testset verbose=true "Wavelength array correctly sorted" begin
        @test issorted(EXOSPHID.solar_wavelength)
    end

    @testset verbose=true "Normalized Solar Flux Database Size" begin
        tol = 1e-6
        for pt in exosphid_species
            fq, fa = get_normalized_fluxes(pt)

            @testset "$pt" begin
                @testset "Array Lengths" begin
                    @test length(fq) == length(EXOSPHID.solar_wavelength)
                    @test length(fa) == length(EXOSPHID.solar_wavelength)
                end

                @testset "Normalization" begin
                    @test abs(sum(fq) - 1.0) < tol
                    @test abs(sum(fa) - 1.0) < tol
                end
            end
        end
    end


    @testset verbose=true "Standard Solar Flux Database Size" begin
        for pt in exosphid_species
            fq, fa = get_standard_fluxes(pt)

            @testset "$pt" begin
                @testset "Array Lengths" begin
                    @test length(fq) == length(EXOSPHID.solar_wavelength)
                    @test length(fa) == length(EXOSPHID.solar_wavelength)
                end
            end
        end
    end

    @testset verbose=true "Normalized Fluxes Compatible with Standard Fluxes" begin

        for pt in exosphid_species
            sq, sa = get_standard_fluxes(pt)
            nq, na = get_normalized_fluxes(pt)

            @testset "$pt" begin
                @test all(isapprox.(nq, EXOSPHID.normalize_flux_distribution(EXOSPHID.solar_wavelength, sq, (1, 95000))[2] ; rtol=1e-6))
                @test all(isapprox.(na, EXOSPHID.normalize_flux_distribution(EXOSPHID.solar_wavelength, sa, (1, 95000))[2] ; rtol=1e-6))
            end
        end
    end

end

println("............... COMPLETED TESTING of solar_database.jl\n")

