using Test
using LinearAlgebra

include("../../src/photoreactions/SimplePhotodissociation.jl")
import .SimplePhotodissociation as SD

include("../../src/database/photodatabase.jl")
using .photodatabase

include("../../src/database/solar_spectrum.jl")
import .solar_spectrum as SS

const eV_to_J = 1.602f-19
const m_fund = 1.66054e-27 # 1 M.U.


@testset verbose=true "SimplePhotodissociation.jl" begin
    @testset verbose=true "random_unit_tuple()" begin
        v = SD.random_unit_tuple()
        
        # Type correctness
        @test v isa Tuple{Float32,Float32,Float32}

        # Magnitude should be 1 (within tolerance)
        @test isapprox(norm(collect(v)), 1f0; atol=1e-5)

        # Randomness (different outputs)
        v2 = SD.random_unit_tuple()
        @test v != v2  # Not guaranteed, but likely true
    end

    @testset verbose=true "calculate_photon_momentum()" begin
        c = 2.99792458f8
        E = 1.0f0
        direction = (1f0, 0f0, 0f0)
        p = SD.calculate_photon_momentum(E, direction)

        # Type check
        @test p isa Tuple{Float32,Float32,Float32}

        # Magnitude check: |p| = E/c
        @test isapprox(norm(collect(p)), E / c; atol=1e-7)
    end


    @testset verbose=true "calculate_excess_energy()" begin
        mp = 17 * m_fund
        Ep = 5.0f0 * eV_to_J
        Eb = 4.0f0 * eV_to_J
        vp = 600.0f0 .* SD.random_unit_tuple()
        st = SD.random_unit_tuple()

        # Case 1: OH but not OH(DPD)
        pn = ("OH(B2Σ)", "O(1S)", "H")
        r1 = SD.PhotoReaction(Eb, vp, st, pn, false)
        E1 = SD.calculate_excess_energy(r1, mp, Ep)
        @test isapprox(E1, Ep - SD.get_electronic_energy_predis("OH(B2Σ)") -SD.get_vibrorotational_energy("OH(B2Σ)"); atol=1e-6)

        # Case 2: OH(DPD)
        pn = ("OH(B2Σ)", "O(1S)", "H")
        r2 = SD.PhotoReaction(Eb, vp, st, pn, false)
        E2 = SD.calculate_excess_energy(r2, mp, Ep)
        @test isapprox(E2, Ep + 0.5f0*mp*norm(vp)^2 - Eb; atol=1e-6)

        # Case 3: Generic reaction (3 products)
        pn = ("H2O", "OH(A2Σ+)", "H")
        r3 = SD.PhotoReaction(Eb, vp, st, pn, false)
        E3 = SD.calculate_excess_energy(r3, mp, Ep)
        E_expected = (Ep + 0.5f0*mp*norm(vp)^2 + SD.get_vibrorotational_energy("H2O")) - Eb - (SD.get_vibrorotational_energy("H") + SD.get_vibrorotational_energy("OH(A2Σ+)"))
        @test isapprox(E3, E_expected; atol=1e-6)
    end

    @testset verbose=true "Test Positive excess energy and velocities" begin

        vps = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, "He" => 1250.0, "Ne" => 560.0)

        num_reactions = 1000
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

                vp = Float32(vps[pt]) .* SD.random_unit_tuple()
                st = SD.random_unit_tuple()

                for (i, rt) in enumerate(rts)

                    if rt in ("SPD", "DPD", "DiPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | Wavelength range: $wr" begin
                                    wvls, energs = SS.flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPD"
                                        tsh = tshs[i] * eV_to_J
                                        pn = sns[i]
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)
                                        p, h, l = r1.product_types
                                        m = get_masses(p, heavy_child_name=h, light_child_name=l)

                                        @testset "Positive Excess Energy" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r1, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Non-zero Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r1, m[1], Ep)
                                                pp = SD.calculate_photon_momentum(Ep, st)
                                                vl, vh = SD.allocate_velocity(r1, Ee, m, pp)
                                            tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
                                            end
                                        end 

                                    elseif rt == "DPD"

                                        # First dissociation

                                        tsh = tshs[i][1] * eV_to_J
                                        pn = sns[i][1]
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)
                                        p, h, l = r1.product_types
                                        m = get_masses(p, heavy_child_name=h, light_child_name=l)

                                        @testset "First Step: Excess Energy" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r1, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        vh_aux = zeros(Float32, 3)

                                        @testset "First Step: Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r1, m[1], Ep)
                                                pp = SD.calculate_photon_momentum(Ep, st)
                                                vl, vh = SD.allocate_velocity(r1, Ee, m, pp)
                                                vh_aux = vh
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
                                            end
                                        end

                                        # Second dissociation
                                        tsh = tshs[i][2] * eV_to_J
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vh_aux)
                                        r2 = SD.PhotoReaction(tsh, vp2, st, pn, false)
                                        p, h, l = r2.product_types
                                        m = get_masses(p, heavy_child_name=h, light_child_name=l)

                                        @testset "Second Step: Excess Energy" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r2, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Second Step: Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r2, m[1], Ep)
                                                pp = SD.calculate_photon_momentum(Ep, st)
                                                vl, vh = SD.allocate_velocity(r2, Ee, m, pp)
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
                                            end
                                        end

                                    elseif rt == "DiPI"
                                        tsh = tshs[i][2] * eV_to_J
                                        pn = map(s -> replace(s, r"\(.*\)" => ""), sns[i][2])
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)
                                        p, h, l = r1.product_types
                                        m = get_masses(p, heavy_child_name=h, light_child_name=l)

                                        @testset "Excess Energy" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r1, m[1], Ep)
                                                @test Ee >= 0
                                            end
                                        end

                                        @testset "Velocity Allocation" begin
                                            for Ep in energs
                                                Ee = SD.calculate_excess_energy(r1, m[1], Ep)
                                                pp = SD.calculate_photon_momentum(Ep, st)
                                                vl, vh = SD.allocate_velocity(r1, Ee, m, pp)
                                                tol = 1e-6
                                                @test (Ee==0 && norm(vh) < tol && norm(vl) < tol) || (Ee != 0 && (norm(vh) !=0  || norm(vl) !=0))
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

        vps = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, "He" => 1250.0, "Ne" => 560.0)

        num_reactions = 1000
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

                vp = Float32(vps[pt]) .* SD.random_unit_tuple()
                st = SD.random_unit_tuple()

                for (i, rt) in enumerate(rts)

                    if rt in ("SPD", "DPD", "DiPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | Wavelength range: $wr" begin
                                    wvls, energs = SS.flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPD"
                                        tsh = tshs[i] * eV_to_J
                                        pn = sns[i]
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)

                                        for (wvl, Ep) in zip(wvls, energs)
                                            if any(r -> r[1] <= wvl <= r[2], wrs)
                                                vl, vh = SD.simulate_photodissociation(r1, Ep) 
                                                @test vl isa NTuple{3, Float32}
                                                @test vh isa NTuple{3, Float32}
                                            end
                                        end
                                        

                                    elseif rt == "DPD"

                                        # First dissociation
                                        tsh = tshs[i][1] * eV_to_J
                                        pn = sns[i][1]
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)
                                        for (wvl, Ep) in zip(wvls, energs)
                                            if any(r -> r[1] <= wvl <= r[2], wrs)
                                                vl, vh = SD.simulate_photodissociation(r1, Ep) 
                                                @test vl isa NTuple{3, Float32}
                                                @test vh isa NTuple{3, Float32}
                                            end
                                        end

                                        # Second dissociation
                                        tsh = tshs[i][2] * eV_to_J
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vh)
                                        r2 = SD.PhotoReaction(tsh, vp2, st, pn, false)
                                        for (wvl, Ep) in zip(wvls, energs)
                                            if any(r -> r[1] <= wvl <= r[2], wrs)
                                                vl, vh = SD.simulate_photodissociation(r2, Ep) 
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

        vps = Dict("H2O" => 590.0, "OH" => 605.0, "H2" => 1750.0, "H" => 2500.0, "H(-)"=> 2500.0, "HO2" => 425.0, "H2O2" => 435.0, "He" => 1250.0, "Ne" => 560.0)

        num_reactions = 1000
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

                vp = Float32(vps[pt]) .* SD.random_unit_tuple()
                st = SD.random_unit_tuple()

                for (i, rt) in enumerate(rts)

                    if rt in ("SPD", "DPD", "DiPI")

                        rn = rns[i]

                        for (wi, wr) in enumerate(wrs)
                            if rps[i][wi] != 0.0
                                @testset verbose=false "Reaction: $rn | Wavelength range: $wr" begin
                                    wvls, energs = SS.flux_outputs(pt, wr, sa, num_reactions)

                                    if rt == "SPD"
                                        tsh = tshs[i] * eV_to_J
                                        pn = sns[i]
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)
                                        vl, vh = SD.multiple_photodissociation(r1, energs)

                                        @test all(x -> x isa NTuple{3, Float32}, vl)
                                        @test length(vl) == num_reactions
                                        @test all(x -> x isa NTuple{3, Float32}, vh)
                                        @test length(vh) == num_reactions

                                    elseif rt == "DPD"

                                        # First dissociation
                                        tsh = tshs[i][1] * eV_to_J
                                        pn = sns[i][1]
                                        r1 = SD.PhotoReaction(tsh, vp, st, pn, false)
                                        vl, vh = SD.multiple_photodissociation(r1, energs)

                                        @test all(x -> x isa NTuple{3, Float32}, vl)
                                        @test length(vl) == num_reactions
                                        @test all(x -> x isa NTuple{3, Float32}, vh)
                                        @test length(vh) == num_reactions

                                        # Second dissociation
                                        tsh = tshs[i][2] * eV_to_J
                                        pn = sns[i][2]
                                        vp2 = map(Float32, vh[end])
                                        r2 = SD.PhotoReaction(tsh, vp2, st, pn, false)
                                        vl, vh = SD.multiple_photodissociation(r2, energs)

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