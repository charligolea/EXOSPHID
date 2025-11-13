using Test
using LinearAlgebra


println("TESTING SimplePhotoionisation.jl ............... ")

@testset verbose=true "SimplePhotoionisation.jl" begin

    @testset verbose=true "random_unit_tuple()" begin
        v = EXOSPHID.random_unit_tuple()
        
        # Type correctness
        @test v isa Tuple{Float32,Float32,Float32}

        # Magnitude should be 1 (within tolerance)
        @test isapprox(norm(collect(v)), 1f0; atol=1e-5)

        # Randomness (different outputs)
        v2 = EXOSPHID.random_unit_tuple()
        @test v != v2  # Not guaranteed, but likely true
    end

    @testset verbose=true "calculate_photon_momentum()" begin
        c = 2.99792458f8
        E = 1.0f0
        direction = (1f0, 0f0, 0f0)
        p = EXOSPHID.calculate_photon_momentum(E, direction)

        # Type check
        @test p isa Tuple{Float32,Float32,Float32}

        # Magnitude check: |p| = E/c
        @test isapprox(norm(collect(p)), E / c; atol=1e-7)
    end

    @testset verbose=true "Test Positive excess energy and velocities" begin

        vps = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, 
                   "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, 
                   "He" => 1250.0, "Ne" => 560.0)

        num_reactions = 10_000
        sa = 0

        for pt in exosphid_species

            @testset verbose=true "$pt" begin
                
                photo_info = get_species_photochemical_info(pt)
                tshs = photo_info.tsh_energies
                sns = photo_info.species_names
                rps = photo_info.reaction_probabilities
                rns = photo_info.reaction_names
                rts = photo_info.reaction_types
                wrs = photo_info.wavelength_range

                vp = Float32(vps[pt]) .* EXOSPHID.random_unit_tuple()
                st = EXOSPHID.random_unit_tuple()
                mp = EXOSPHID.get_masses(pt)

                for (i, rt) in enumerate(rts)

                    if rt in ("SPI", "DPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | λ: $wr" begin
                                    wvls, energs = flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPI"
                                        tsh = EXOSPHID.eV2J(tshs[i])
                                        pn = sns[i]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)

                                        @testset "Positive Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_ionisation(tsh, Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Non-zero Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_ionisation(tsh, Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vion = EXOSPHID.allocate_velocity_ionisation(r1, Ee, mp, pp)
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vion)) || 
                                                      (Ee != 0 && (norm(vion) !=0))
                                            end
                                        end 

                                    elseif rt == "DPI"

                                        # First dissociation

                                        tsh = EXOSPHID.eV2J(tshs[i][1])
                                        pn = sns[i][1]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)

                                        @testset "First Step: Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_ionisation(tsh, Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        vion_aux = zeros(Float32, 3)

                                        @testset "First Step: Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_ionisation(tsh, Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vion = EXOSPHID.allocate_velocity_ionisation(r1, Ee, mp, pp)
                                                vion_aux = vion
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vion) < tol) || 
                                                      (Ee != 0 && (norm(vion) !=0))
                                            end
                                        end

                                        # Second dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vion_aux)
                                        r2 = EXOSPHID.PhotoReaction(tsh, vp2, st, pn, false)

                                        @testset "Second Step: Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_ionisation(tsh, Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Second Step: Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_ionisation(tsh, Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vion = EXOSPHID.allocate_velocity_ionisation(r2, Ee, mp, pp)
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vion) < tol) || 
                                                    (Ee != 0 && norm(vion) !=0)
                                            end
                                        end

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    @testset verbose=true "simulate_photoionisation() output type" begin

        vps = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, 
                   "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, 
                   "He" => 1250.0, "Ne" => 560.0)

        num_reactions = 10_000
        sa = 0

        for pt in exosphid_species

            @testset verbose=true "$pt" begin
                
                photo_info = get_species_photochemical_info(pt)
                tshs = photo_info.tsh_energies
                sns = photo_info.species_names
                rps = photo_info.reaction_probabilities
                rns = photo_info.reaction_names
                rts = photo_info.reaction_types
                wrs = photo_info.wavelength_range

                vp = Float32(vps[pt]) .* EXOSPHID.random_unit_tuple()
                st = EXOSPHID.random_unit_tuple()
                mp = EXOSPHID.get_masses(pt)

                for (i, rt) in enumerate(rts)

                    if rt in ("SPI", "DPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | λ: $wr" begin
                                    wvls, energs = flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPI"
                                        tsh = EXOSPHID.eV2J(tshs[i])
                                        pn = sns[i]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)

                                        for Ep in energs
                                            @test simulate_photoionisation(r1, Ep) isa 
                                                NTuple{3, Float32}
                                        end

                                    elseif rt == "DPI"

                                        # First dissociation

                                        tsh = EXOSPHID.eV2J(tshs[i][1])
                                        pn = sns[i][1]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)

                                        vion = zeros(Float32, 3)

                                        for Ep in energs
                                            vion = simulate_photoionisation(r1, Ep)
                                            @test vion isa NTuple{3, Float32}
                                        end


                                        # Second dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vion)
                                        r2 = EXOSPHID.PhotoReaction(tsh, vp2, st, pn, false)

                                        for Ep in energs
                                            @test simulate_photoionisation(r2, Ep) isa 
                                                    NTuple{3, Float32}
                                        end

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end

    @testset verbose=true "multiple_photoionisation() output type and size" begin

        vps = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, 
                   "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, 
                   "He" => 1250.0, "Ne" => 560.0)

        num_reactions = 10_000
        sa = 0

        for pt in exosphid_species

            @testset verbose=true "$pt" begin
                
                photo_info = get_species_photochemical_info(pt)
                tshs = photo_info.tsh_energies
                sns = photo_info.species_names
                rps = photo_info.reaction_probabilities
                rns = photo_info.reaction_names
                rts = photo_info.reaction_types
                wrs = photo_info.wavelength_range

                vp = Float32(vps[pt]) .* EXOSPHID.random_unit_tuple()
                st = EXOSPHID.random_unit_tuple()
                mp = EXOSPHID.get_masses(pt)

                for (i, rt) in enumerate(rts)

                    if rt in ("SPI", "DPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | : $wr" begin
                                    wvls, energs = flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPI"
                                        tsh = EXOSPHID.eV2J(tshs[i])
                                        pn = sns[i]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        vel = multiple_photoionisation(r1, energs) 
                                        @test all(x -> x isa NTuple{3, Float32}, vel)
                                        @test length(vel) == num_reactions

                                    elseif rt == "DPI"

                                        # First dissociation

                                        tsh = EXOSPHID.eV2J(tshs[i][1])
                                        pn = sns[i][1]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)

                                        vel = multiple_photoionisation(r1, energs) 
                                        @test all(x -> x isa NTuple{3, Float32}, vel)
                                        @test length(vel) == num_reactions


                                        # Second dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vel[end])
                                        r2 = EXOSPHID.PhotoReaction(tsh, vp2, st, pn, false)

                                        vel = multiple_photoionisation(r2, energs) 
                                        @test all(x -> x isa NTuple{3, Float32}, vel)
                                        @test length(vel) == num_reactions

                                    end

                                end

                            end

                        end

                    end

                end

            end

        end

    end


end

println("............... COMPLETED TESTING of SimplePhotoionisation.jl\n")
