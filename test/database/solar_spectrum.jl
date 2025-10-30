println("TESTING solar_spectrum.jl ...............")

@testset verbose=true "solar_spectrum.jl" begin

    # Even if you select a smaller range of wvl check that sorted and same size

    for pt in exosphid_species

        @testset "$pt" verbose=true begin 

            for i in 1:5 # Try out 1000 cases

                wvr1, wvr2 = sort(rand(1:95_000, 2))
                tol = 1e-5
                
                wvl, fq, fa = get_solar_fluxes(pt, (wvr1, wvr2))
                
                @testset "Check flux Array Lengths are correct for reduced WVL range" begin
                    @test length(fq) == length(wvl)
                    @test length(fa) == length(wvl)
                end

                @testset "Check Normalization is well done for reduced WVL range" begin
                    @test abs(sum(fq) - 1.0) < tol
                    @test abs(sum(fa) - 1.0) < tol
                end 
            end
        end
    end
end


println("............... COMPLETED TESTING solar_spectrum.jl\n")