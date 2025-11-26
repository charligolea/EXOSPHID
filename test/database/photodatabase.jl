println("TESTING photodatabase.jl ...............")


# ─────────────────────────────────────────────────────────────────────────────────────────
# TEST 1: Validate Structure of Photochemical Database
# ─────────────────────────────────────────────────────────────────────────────────────────

@testset verbose=true "Photodatabase Tables Validation" begin

    # 0. Test that database exists
    @testset verbose=true "Database exists for all parents in exosphid_species" begin
        @test all(x -> isfile(joinpath(@__DIR__, "..", "..", "src", "database", "species", "$x.jl")), exosphid_species)
    end

    for parent_type in exosphid_species

        @testset verbose=false "$parent_type database" begin

            photo_info = get_species_photochemical_info(parent_type)
            
            # 1. Test that none of the returned tables are nothing
            @testset verbose=false "Database is complete" begin
                @test !isnothing(photo_info.tsh_energies)
                @test !isnothing(photo_info.species_names)
                @test !isnothing(photo_info.reaction_probabilities)
                @test !isnothing(photo_info.reaction_names)
                @test !isnothing(photo_info.reaction_types)
                @test !isnothing(photo_info.wavelength_range)
            end
            
            N = length(photo_info.reaction_types)

            # 2. All tables have the right length
            @testset verbose=false "Test data has the right length" begin
                @test length(photo_info.tsh_energies) == N
                @test length(photo_info.species_names) == N
                @test length(photo_info.reaction_probabilities) == N
                @test length(photo_info.reaction_names) == N
            end

            # 3. Each element has correct information, type and structure
            @testset verbose=false "Database has correct structure and type" begin
                for (i, rt) in enumerate(photo_info.reaction_types)
                    
                    rn = photo_info.reaction_names[i]
                    @testset verbose=false "$rn" begin

                        @test rt in ("SPD", "SPI", "DPD", "DPI", "DiPI") 

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

                                sns = map(s -> replace(s, r"\(.*\)" => ""),sn)
                                if rt == "SPI"
                                    @test sn[1] in exosphid_species
                                    if !occursin("(-)", sn[1]) # negative hydrogen case
                                        @test sn[2] == sn[1] * "(+)"
                                    else
                                        @test sn[2] == sns[1]
                                    end
                                    @test sn[3] == "e(-)"
                                else
                                    masses = (EXOSPHID.get_masses(sns[1]), 
                                              EXOSPHID.get_masses(sns[2]), 
                                              EXOSPHID.get_masses(sns[3]))
                                    # This checks that the species are correctly ordered: 
                                    # parent - heavy - light
                                    # It may be that there is no heavy and light product, 
                                    # such as H2 -> H + H
                                    @test masses[1] > masses[2] && masses[2]  >= masses[3]
                                end

                            elseif rt in ("DPD", "DPI", "DiPI")
                                @test sn isa Tuple
                                @test length(sn) == 2
                                @test all(sub -> sub isa Tuple && length(sub) == 3 && 
                                      all(x -> x isa String, sub), sn)
                                
                                if rt == "DPI"

                                    for species_tuple in sn
                                        @test species_tuple[1] in exosphid_species
                                        # Same checks as SPI!!!
                                        sns = map(s -> replace(s, r"\(.*\)" => ""),
                                                  species_tuple)
                                        
                                        # negative hydrogen case
                                        if !occursin("(-)", species_tuple[1]) 
                                            @test species_tuple[2] == 
                                                    species_tuple[1] * "(+)"
                                        else
                                            @test species_tuple[2] == sns[1]
                                        end
                                        @test species_tuple[3] == "e(-)"
                                    end

                                elseif rt == "DPD"
                                    for st in sn
                                        species_tuple = 
                                            map(s -> replace(s, r"\(.*\)" => ""),st)
                                        masses = (EXOSPHID.get_masses(species_tuple[1]), 
                                                  EXOSPHID.get_masses(species_tuple[2]), 
                                                  EXOSPHID.get_masses(species_tuple[3]))
                                        @test masses[1] > masses[2] && 
                                              masses[2]  >= masses[3]
                                    end

                                elseif rt == "DiPI"
                                    # Combine tests above

                                    # Ionisation tests for sn[1]
                                    @test sn[1][1] in exosphid_species
                                    if !occursin("(-)", sn[1][1]) # negative hydrogen case
                                        @test sn[1][2] == sn[1][1] * "(+)"
                                    else
                                        @test sn[1][2] == sns[1][1]
                                    end
                                    @test sn[1][3] == "e(-)"

                                    # Dissociation tests for sn[2]
                                    sns = map(s -> replace(s, r"\(.*\)" => ""),sn[2])
                                    masses = (EXOSPHID.get_masses(sns[1]), 
                                              EXOSPHID.get_masses(sns[2]), 
                                              EXOSPHID.get_masses(sns[3]))
                                    @test masses[1] > masses[2] && masses[2]  >= masses[3]
                                end

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
                @test all(range -> range[2] <= 
                            photo_info.wvl_threshold, photo_info.wavelength_range)
            end

            # 6. Reaction probabilities for specific wavelength range sum 1
            @testset verbose=true "Cumulated probability per wvl range is 1" begin
                @test all(sum(photo_info.reaction_probabilities[j][i] for j in 1:N) == 1.0 
                            for i in 1:length(photo_info.wavelength_range))
            end

            # 7. Positive photodestruction rates
            @testset verbose=true "Postive Photodestruction rates" begin
                @test 0 < photo_info.quiet_rate <= photo_info.active_rate
                @test EXOSPHID.get_photodestruction_rates(photo_info, rand(Float32), 1.0) > 0
            end

            # 8. Test that active photo rates are larger or equal to quiet photo rates
            @testset verbose=true "Active Photodestruction rate >= quiet rate" begin
                @test photo_info.quiet_rate <= photo_info.active_rate
            end

            # 9. Make sure photo rates are consistent to the solar flux database
            @testset verbose=false "Photorates consistent with fluxes" begin

                @testset verbose=false "Normalized" begin
                    qf, af = get_normalized_fluxes(parent_type)

                    if photo_info.quiet_rate != photo_info.active_rate
                        @test any(.!isapprox.(qf, af; rtol=1e-12))
                    elseif photo_info.quiet_rate == photo_info.active_rate
                        all(isapprox.(qf, af; rtol=1e-12))
                    end
                end

                @testset verbose=false "Standard" begin
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


# ─────────────────────────────────────────────────────────────────────────────────────────
# TEST 2: Test various functions from photodatabase.jl
# ─────────────────────────────────────────────────────────────────────────────────────────

@testset verbose = true "is_photoreaction_occuring" begin
    @test EXOSPHID.is_photoreaction_occuring(rand(Float32), rand(Float32)) isa Bool
end

@testset verbose = true "get_current_reaction" begin

    for parent_type in exosphid_species
        @testset verbose=false "$parent_type database" begin

            photo_info = get_species_photochemical_info(parent_type)

            @testset verbose=false "Photons with wvl in relevant range trigger reaction" begin
            
                p_wvl = photo_info.wvl_threshold * (1.f0+rand(Float32)) # Photon wavelength higher than threshold
                @testset verbose=false "No reaction if wvl > threshold" begin
                    @test get_current_reaction(p_wvl, photo_info) == ("", "", Integer[], ())
                end

                p_wvl = rand() * photo_info.wvl_threshold # Photon wavelength lower than threshold
                @testset verbose=false "Yes reaction if wvl < threshold AND in relevant wvl_range" begin
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

@testset verbose = false "Species masses" begin
    @testset verbose=false "Mass available for given species" begin
        @test length(EXOSPHID.mass_species) == length(EXOSPHID.mass_dict)
    end
    @testset verbose=false "Masses are positive and within expected range" begin
        @test all(0 .<= EXOSPHID.mass_dict .<= EXOSPHID.amu2kg(100))
    end
    @testset "get_masses returns numeric values" begin
        for sp in exosphid_species
            masses = EXOSPHID.get_masses(sp)
            @test masses isa Float64
        end
    end
end


# ─────────────────────────────────────────────────────────────────────────────────────────
# TEST 3: Test vibrorotational and electronic energy database
# ─────────────────────────────────────────────────────────────────────────────────────────

@testset verbose = false "Vibrorotational energy" begin

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
                            @testset verbose=false "$sp" begin
                                @testset verbose=false "Available" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) !== nothing
                                end
                                @testset verbose=false "Positive" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) >=0.0
                                end
                            end
                        end
                    end
                elseif rt == "DPD"
                    for sp in spns[i][1]
                        if sp ∉ tested_states
                            push!(tested_states, sp)
                            @testset verbose=false "$sp" begin
                                @testset verbose=false "Available" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) !== nothing
                                end
                                @testset verbose=false "Positive" begin
                                    @test EXOSPHID.get_vibrorotational_energy(sp) >= 0.0
                                end
                            end
                        end
                    end
                elseif rt == "DiPI"
                    for sp in map(s -> replace(s, r"\(.*\)" => ""),(spns[i][2]))
                        if sp ∉ tested_states
                            push!(tested_states, sp)
                            @testset verbose=false "$sp" begin
                                @testset verbose=false "Available" begin 
                                    @test EXOSPHID.get_vibrorotational_energy(sp) !== nothing
                                end
                                @testset verbose=false "Positive" begin 
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

@testset verbose = false "Electronic predissociation energy" begin

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
                    @testset verbose=false "$sp_OH" begin
                        @testset verbose=false "Available" begin
                            @test EXOSPHID.get_electronic_energy_predis(sp_OH) !== nothing
                        end
                        @testset verbose=false "Positive" begin
                            @test EXOSPHID.get_electronic_energy_predis(sp_OH) >=0.0
                        end
                    end
                end
            end
        end
    end
end


println("............... COMPLETED TESTING photodatabase.jl\n")