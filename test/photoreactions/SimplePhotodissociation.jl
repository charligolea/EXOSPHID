using Test
using LinearAlgebra

const eV_to_J = 1.602f-19
const m_fund = 1.66054e-27 # 1 M.U.

println("TESTING SimplePhotodissociation.jl ............... ")

@testset verbose=true "SimplePhotodissociation.jl" begin
    @testset verbose=true "EXOSPHID.random_unit_tuple()" begin
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


    @testset verbose=true "calculate_excess_energy_dissociation()" begin
        mp = 17 * m_fund
        Ep = EXOSPHID.eV2J(5.0f0)
        Eb = EXOSPHID.eV2J(4.0f0)
        vp = 600.0f0 .* EXOSPHID.random_unit_tuple()
        st = EXOSPHID.random_unit_tuple()

        # Case 1: OH but not OH(DPD)
        pn = ("OH(B2Σ)", "O(1S)", "H")
        r1 = EXOSPHID.PhotoReaction(Eb, vp, st, pn, false)
        E1 = EXOSPHID.calculate_excess_energy_dissociation(r1, mp, Ep)
        @test isapprox(E1, Ep - EXOSPHID.get_electronic_energy_predis("OH(B2Σ)") -
                        EXOSPHID.get_vibrorotational_energy("OH(B2Σ)"); atol=1e-6)

        # Case 2: OH(DPD)
        pn = ("OH(B2Σ)", "O(1S)", "H")
        r2 = EXOSPHID.PhotoReaction(Eb, vp, st, pn, false)
        E2 = EXOSPHID.calculate_excess_energy_dissociation(r2, mp, Ep)
        @test isapprox(E2, Ep + 0.5f0*mp*norm(vp)^2 - Eb; atol=1e-6)

        # Case 3: Generic reaction (3 products)
        pn = ("H2O", "OH(A2Σ+)", "H")
        r3 = EXOSPHID.PhotoReaction(Eb, vp, st, pn, false)
        E3 = EXOSPHID.calculate_excess_energy_dissociation(r3, mp, Ep)
        E_expected = (Ep + 0.5f0*mp*norm(vp)^2 + EXOSPHID.get_vibrorotational_energy("H2O")) 
                        - Eb - (EXOSPHID.get_vibrorotational_energy("H") + 
                                EXOSPHID.get_vibrorotational_energy("OH(A2Σ+)"))
        @test isapprox(E3, E_expected; atol=1e-6)
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

                for (i, rt) in enumerate(rts)

                    if rt in ("SPD", "DPD", "DiPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | λ: $wr" begin
                                    wvls, energs = flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPD"
                                        tsh = EXOSPHID.eV2J(tshs[i])
                                        pn = sns[i]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        p, h, l = r1.product_types
                                        m = (EXOSPHID.get_masses(p), 
                                             EXOSPHID.get_masses(h), 
                                             EXOSPHID.get_masses(l))

                                        @testset "Positive Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r1, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Non-zero Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r1, m[1], Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vl, vh = EXOSPHID.allocate_velocity_dissociation(r1, Ee, m, pp)
                                            tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || 
                                                      (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
                                            end
                                        end 

                                    elseif rt == "DPD"

                                        # First dissociation

                                        tsh = EXOSPHID.eV2J(tshs[i][1])
                                        pn = sns[i][1]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        p, h, l = r1.product_types
                                        m = (EXOSPHID.get_masses(p), 
                                             EXOSPHID.get_masses(h), 
                                             EXOSPHID.get_masses(l))

                                        @testset "First Step: Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r1, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        vh_aux = zeros(Float32, 3)

                                        @testset "First Step: Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r1, m[1], Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vl, vh = EXOSPHID.allocate_velocity_dissociation(r1, Ee, m, pp)
                                                vh_aux = vh
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || 
                                                      (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
                                            end
                                        end

                                        # Second dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vh_aux)
                                        r2 = EXOSPHID.PhotoReaction(tsh, vp2, st, pn, false)
                                        p, h, l = r2.product_types
                                        m = (EXOSPHID.get_masses(p), 
                                             EXOSPHID.get_masses(h), 
                                             EXOSPHID.get_masses(l))

                                        @testset "Second Step: Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r2, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Second Step: Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r2, m[1], Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vl, vh = EXOSPHID.allocate_velocity_dissociation(r2, Ee, m, pp)
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || 
                                                      (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
                                            end
                                        end

                                    elseif rt == "DiPI"
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = map(s -> replace(s, r"\(.*\)" => ""), sns[i][2])
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        p, h, l = r1.product_types
                                        m = (EXOSPHID.get_masses(p), 
                                             EXOSPHID.get_masses(h), 
                                             EXOSPHID.get_masses(l))

                                        @testset "Excess Energy" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r1, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = EXOSPHID.calculate_excess_energy_dissociation(r1, m[1], Ep)
                                                pp = EXOSPHID.calculate_photon_momentum(Ep, st)
                                                vl, vh = EXOSPHID.allocate_velocity_dissociation(r1, Ee, m, pp)
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || 
                                                    (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
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

    @testset verbose=true "simulate_photodissociation() output type" begin

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

                for (i, rt) in enumerate(rts)

                    if rt in ("SPD", "DPD", "DiPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | λ: $wr" begin
                                    wvls, energs = flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPD"
                                        tsh = EXOSPHID.eV2J(tshs[i])
                                        pn = sns[i]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)

                                        for (wvl, Ep) in zip(wvls, energs)
                                            if any(r -> r[1] <= wvl <= r[2], wrs)
                                                vl, vh = simulate_photodissociation(r1, Ep) 
                                                @test vl isa NTuple{3, Float32}
                                                @test vh isa NTuple{3, Float32}
                                            end
                                        end
                                        

                                    elseif rt == "DPD"

                                        vh = nothing 

                                        # First dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][1])
                                        pn = sns[i][1]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        for (wvl, Ep) in zip(wvls, energs)
                                            if any(r -> r[1] <= wvl <= r[2], wrs)
                                                vl, vh = simulate_photodissociation(r1, Ep) 
                                                @test vl isa NTuple{3, Float32}
                                                @test vh isa NTuple{3, Float32}
                                            end
                                        end

                                        # Second dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vh)
                                        r2 = EXOSPHID.PhotoReaction(tsh, vp2, st, pn, false)
                                        for (wvl, Ep) in zip(wvls, energs)
                                            if any(r -> r[1] <= wvl <= r[2], wrs)
                                                vl, vh = simulate_photodissociation(r2, Ep) 
                                                @test vl isa NTuple{3, Float32}
                                                @test vh isa NTuple{3, Float32}
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

    @testset verbose=true "multiple_photodissociation() output type and length" begin

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

                for (i, rt) in enumerate(rts)

                    if rt in ("SPD", "DPD", "DiPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | λ: $wr" begin
                                    wvls, energs = flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPD"
                                        tsh = EXOSPHID.eV2J(tshs[i])
                                        pn = sns[i]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        vl, vh = multiple_photodissociation(r1, energs)

                                        @test all(x -> x isa NTuple{3, Float32}, vl)
                                        @test length(vl) == num_reactions
                                        @test all(x -> x isa NTuple{3, Float32}, vh)
                                        @test length(vh) == num_reactions

                                    elseif rt == "DPD"

                                        # First dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][1])
                                        pn = sns[i][1]
                                        r1 = EXOSPHID.PhotoReaction(tsh, vp, st, pn, false)
                                        vl, vh = multiple_photodissociation(r1, energs)

                                        @test all(x -> x isa NTuple{3, Float32}, vl)
                                        @test length(vl) == num_reactions
                                        @test all(x -> x isa NTuple{3, Float32}, vh)
                                        @test length(vh) == num_reactions

                                        # Second dissociation
                                        tsh = EXOSPHID.eV2J(tshs[i][2])
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vh[end])
                                        r2 = EXOSPHID.PhotoReaction(tsh, vp2, st, pn, false)
                                        vl, vh = multiple_photodissociation(r2, energs)

                                        @test all(x -> x isa NTuple{3, Float32}, vl)
                                        @test length(vl) == num_reactions
                                        @test all(x -> x isa NTuple{3, Float32}, vh)
                                        @test length(vh) == num_reactions

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

println("............... COMPLETED TESTING of SimplePhotodissociation.jl\n")
