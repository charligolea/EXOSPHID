const m_fund = 1.66054e-27 # 1 M.U.

println("TESTING photodatabase.jl ...............")

@testset verbose=true "Photodatabase Tables Validation" begin

    # 0. Test that database exists
    @testset verbose=true "Database exists for all parents in exosphid_species" begin
        for parent in exosphid_species
            @test isfile(joinpath(@__DIR__, "..", "..", "src", "database", "species", "$parent.jl"))
        end
    end

    for parent_type in exosphid_species

        @testset verbose=true "$parent_type database" begin

            photo_info = get_species_photochemical_info(parent_type)
            
            # 1. Test that none of the returned tables are nothing
            @testset verbose=true "Database is complete" begin
                @test photo_info.tsh_energies !== nothing
                @test photo_info.species_names !== nothing
                @test photo_info.reaction_probabilities !== nothing
                @test photo_info.reaction_names !== nothing
                @test photo_info.reaction_types !== nothing
                @test photo_info.wavelength_range !== nothing
            end
            
            N = length(photo_info.reaction_types)

            # 2. All tables have the same length
            @testset verbose=true "Test data has the right length" begin
                @test length(photo_info.tsh_energies) == N
                @test length(photo_info.species_names) == N
                @test length(photo_info.reaction_probabilities) == N
                @test length(photo_info.reaction_names) == N
            end

            # 3. Each element
            @testset verbose=true "Database has correct structure and type" begin
                for (i, rt) in enumerate(photo_info.reaction_types)
                    
                    rn = photo_info.reaction_names[i]
                    @testset verbose=true "$rn" begin

                        @testset verbose=true "Acceptable reaction type" begin
                            @test rt in ("SPD", "SPI", "DPD", "DPI", "DiPI") 
                        end

                        # a) tsh_energies
                        @testset verbose=true "Threshold energies" begin
                            te = photo_info.tsh_energies[i]
                            if rt in ("SPD", "SPI")
                                @test te isa Float32 || te isa Float64
                            elseif rt in ("DPD", "DPI", "DiPI")
                                @test te isa Tuple
                                @test length(te) == 2
                                @test all(x -> x isa Float32 || x isa Float64, te)
                                @test te[2]>te[1]
                            end
                        end

                        # b) species_names structure
                        @testset verbose=true "Species Names" begin
                            sn = photo_info.species_names[i]
                            if rt in ("SPD", "SPI")
                                @test sn isa Tuple
                                @test length(sn) == 3
                                @test all(x -> x isa String, sn)
                            elseif rt in ("DPD", "DPI", "DiPI")
                                @test sn isa Tuple
                                @test length(sn) == 2
                                @test all(sub -> sub isa Tuple && length(sub) == 3 && all(x -> x isa String, sub), sn)
                            end
                        end

                        # c) reaction_probabilities length
                        @testset verbose=true "Reaction probabilities" begin
                            rp = photo_info.reaction_probabilities[i]
                            @test length(rp) == length(photo_info.wavelength_range)
                        end
                    
                    end
                end

                @testset verbose=true "Wavelength range" begin
                    @test all(length.(photo_info.wavelength_range) .== 2)
                    @test all(r -> r[2] > r[1], photo_info.wavelength_range)
                end

            end

            # 4. Wavelength threshold
            @testset verbose=true "Wavelength threshold is positive" begin
                @test photo_info.wvl_threshold isa Float32
                @test photo_info.wvl_threshold > 0
            end

            # 5. Species spefic Wavelength Ranges are always below threshold
            @testset verbose=true "Wavelength ranges are below threshold" begin
                @test all(range -> range[2] <= photo_info.wvl_threshold, photo_info.wavelength_range)
            end

            # 6. Reaction probabilities for specific wavelength range sum 1
            @testset verbose=true "Cumulated probability per wvl range is 1" begin
                @test all(sum(photo_info.reaction_probabilities[j][i] for j in 1:N) == 1.0 for i in 1:length(photo_info.wavelength_range))
            end

            # 7. Positive photodestruction rates
            @testset verbose=true "Postive Photodestruction rates" begin
                @test photo_info.quiet_rate > 0
                @test photo_info.quiet_rate > 0
                @test EXOSPHID.get_photodestruction_rates(photo_info, rand(Float32), 1.0) > 0
            end

            @testset verbose=true "Photorates consistent with fluxes" begin

                @testset verbose=true "Normalized" begin
                    qf, af = get_normalized_fluxes(parent_type)

                    if photo_info.quiet_rate != photo_info.active_rate
                        @test any(.!isapprox.(qf, af; rtol=1e-12))
                    elseif photo_info.quiet_rate == photo_info.active_rate
                        all(isapprox.(qf, af; rtol=1e-12))
                    end
                end

                @testset verbose=true "Standard" begin
                    qf, af = get_standard_fluxes(parent_type)

                    if photo_info.quiet_rate != photo_info.active_rate
                        @test any(.!isapprox.(qf, af; rtol=1e-12))
                    elseif photo_info.quiet_rate == photo_info.active_rate
                        all(isapprox.(qf, af; rtol=1e-12))
                    end
                end

            end

        end
    end
end


@testset verbose = true "is_photoreaction_occuring" begin
    @test EXOSPHID.is_photoreaction_occuring(rand(Float32), rand(Float32)) isa Bool
end

@testset verbose = true "get_current_reaction" begin

    for parent_type in exosphid_species
        @testset verbose=true "$parent_type database" begin

            photo_info = get_species_photochemical_info(parent_type)

            @testset verbose=true "Photons with wvl in relevant range trigger reaction" begin
            
                p_wvl = photo_info.wvl_threshold * (1.f0+rand(Float32)) # Photon wavelength higher than threshold
                @testset verbose=true "No reaction if wvl > threshold" begin
                    @test get_current_reaction(p_wvl, photo_info) == ("", "", Integer[], ())
                end

                p_wvl = rand() * photo_info.wvl_threshold # Photon wavelength lower than threshold
                @testset verbose=true "Yes reaction if wvl < threshold AND in relevant wvl_range" begin
                    if any(r -> r[1] <= p_wvl <= r[2], photo_info.wavelength_range)
                        @test get_current_reaction(p_wvl, photo_info) != ("", "", Integer[], ())
                    else
                        @test get_current_reaction(p_wvl, photo_info) == ("", "", Integer[], ())
                    end
                end
            end
        end
    end
    
end

@testset verbose = true "Vibrorotational energy" begin

    tested_states = []

    for parent_type in exosphid_species
        photo_info = get_species_photochemical_info(parent_type) 
        spns = photo_info.species_names
        rts = photo_info.reaction_types

        for (i,rt) in enumerate(rts)
            if rt in ("SPD", "DPD", "DiPI") 
                if rt == "SPD"
                    for sp in spns[i]
                        if sp ∉ tested_states
                            push!(tested_states, sp)
                            @testset verbose=true "$sp" begin
                                @testset verbose=true "Available" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) !== nothing
                                end
                                @testset verbose=true "Positive" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) >=0.0
                                end
                            end
                        end
                    end
                elseif rt == "DPD"
                    for sp in spns[i][1]
                        if sp ∉ tested_states
                            push!(tested_states, sp)
                            @testset verbose=true "$sp" begin
                                @testset verbose=true "Available" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) !== nothing
                                end
                                @testset verbose=true "Positive" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) >= 0.0
                                end
                            end
                        end
                    end
                elseif rt == "DiPI"
                    for sp in map(s -> replace(s, r"\(.*\)" => ""),(spns[i][2]))
                        if sp ∉ tested_states
                            push!(tested_states, sp)
                            @testset verbose=true "$sp" begin
                                @testset verbose=true "Available" begin 
                                    @test EXOSPHID.get_vibrorotational_energy(sp) !== nothing
                                end
                                @testset verbose=true "Positive" begin 
                                    @test EXOSPHID.get_vibrorotational_energy(sp) >= 0.0
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset verbose = true "Electronic predissociation energy" begin

    tested_states = []

    for parent_type in ("OH", ) # Modify if more species with predissociation are added
        photo_info = get_species_photochemical_info(parent_type) 
        spns = photo_info.species_names
        rts = photo_info.reaction_types

        for (i,rt) in enumerate(rts)
            if rt == "SPD"
                sp_OH = spns[i][1]
                if sp_OH ∉ tested_states
                    push!(tested_states, sp_OH)
                    @testset verbose=true "$sp_OH" begin
                        @testset verbose=true "Available" begin
                            @test EXOSPHID.get_electronic_energy_predis(sp_OH) !== nothing
                        end
                        @testset verbose=true "Positive" begin
                            @test EXOSPHID.get_electronic_energy_predis(sp_OH) >=0.0
                        end
                    end
                end
            end
        end
    end
end


@testset verbose = true "SPECIES MASSES" begin
    @testset verbose=true "Mass available for given species" begin
        @test length(EXOSPHID.mass_species) == length(EXOSPHID.mass_dict)
    end
    @testset verbose=true "Masses are positive and within expected range" begin
        @test all(0 .<= EXOSPHID.mass_dict .<= 300*EXOSPHID.m_fund)
    end
    @testset "get_masses returns numeric values" begin
        for sp in exosphid_species
            masses = EXOSPHID.get_masses(sp, mode="PI")
            @test masses !== nothing
            @test masses isa Float64
        end
    end
end

println("............... COMPLETED TESTING photodatabase.jl\n")