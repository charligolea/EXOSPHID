module ExosPhID

using Random, LinearAlgebra, Distributions

include("photoreactions/SimplePhotodissociation.jl")
using .SimplePhotodissociation

include("photoreactions/SimplePhotoionisation.jl")
using .SimplePhotoionisation

include("solar_spectrum.jl")
using .solar_spectrum

export photodestruction

function get_photodestruction_rates_and_lifetimes(species::String, solar_activity::Float64; mode::String="rate")
    """
    Return the photodestruction rate `k` [1/s] or lifetime `1/k` [s] 
    for a given `species` and `solar_activity` level (0 = quiet, 1 = active).

    Valid species values: "H2O", "OH", "H2", "H", "H(-)"
    Modes: "rate" (default), "lifetime", or "both"
    """

    quiet_rate = 0.0
    active_rate = 0.0

    if species == "H2O"
        quiet_rate = 12.94e-6   # Gomez de Olea (2026)
        active_rate = 22.25e-6
    elseif species == "OH"
        quiet_rate = 22.45e-6   # Gomez de Olea (2025)
        active_rate = 32.56e-6
    elseif species == "H2"
        quiet_rate = 14.24e-6   # Gomez de Olea (2025)
        active_rate = 32.83e-6
    elseif species == "H"
        quiet_rate = 7.00e-8    # Gomez de Olea (2026)
        active_rate = 1.68e-7
    elseif species == "H(-)"
        quiet_rate = 14.24      # Gomez de Olea (2026)
        active_rate = 14.24
    elseif species == "HO2"
        quiet_rate = 6.59-3    # Gomez de Olea (2026)
        active_rate = 6.59-3
    elseif species == "H2O2"
        quiet_rate = 1.43e-4    # Gomez de Olea (2026)
        active_rate = 1.43e-4
    elseif species == "He"
        quiet_rate = 5.57e-8    # Gomez de Olea (2026)
        active_rate = 16.80e-8
    elseif species == "Ne"
        quiet_rate = 1.72e-7    # Gomez de Olea (2026)
        active_rate = 5.74e-7
    else
        throw(ArgumentError("Unknown species: $species"))
    end

    k = (1 - solar_activity) * quiet_rate + solar_activity * active_rate

    return mode == "rate" ? k :
           mode == "lifetime" ? 1/k :
           mode == "both" ? (rate=k, lifetime=1/k) :
           throw(ArgumentError("Mode must be 'rate', 'lifetime' or 'both'"))
end


function photodestruction(solar_activity::Float64, dt::Union{Float64, Int64}, parent_type::String, parent_velocity::Union{Tuple{Float64, Float64, Float64}, Float64}, sun_tuple::Union{Tuple{Float64, Float64, Float64}, Nothing})

    """ 
    INPUTS: 
    - solar_activity = Float from 0.0 (Quiet Sun) to 1.0 (Active Sun)
    - dt: analyzed simulation time window
    - parent_type:: String with parent molecule type ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2")
    - parent_velocity:: Tuple{Float64, Float64, Float64} or Float64, Parent velocity vector in m/s    
    
    OUTPUTS:
    - reaction_occurence::Boolean
    - present_reaction:: String poiting at the reaction that is occuring
    - product_types:: ["T1", "T2", "T3"]
    - product_velocities:: [(vx1, vy1, vz1), (vx2, vy2, vz2), (vx3, vy3, vz3)]
    """

    # Threshold energies are taken from Huebner (2015):
    wvl_thresholds = (1860.0, 2460.0, 2768.0, 911.75, 16439.0, 4395.0, 5765.0, 504.27, 574.94) # "H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne" in Armstrong

    # Initialise variables
    product_types = []
    product_velocities = []
    present_reaction = ""
    reaction_name = ""
    wvl_range = ()

    # 0. Check if reaction is happening according to probability

    k = get_photodestruction_rates_and_lifetimes(parent_type, solar_activity, mode="rate") # Photodestruction rate
    
    p = 1-exp(-k*dt)
    reaction_occurence = rand() < p

    if reaction_occurence == true

        # 1. Generate Photon Energy (J) and Wavelength (A) from Huebner Solar Flux Distributions (1992, 2015)
        photon_energy, photon_wvl = flux_outputs(parent_type, (1.0, 95000.0), solar_activity, 1)

        # 2. Simulate photoreaction

        if parent_type in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne")
            type_index = findfirst(x -> x == parent_type, ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne"))

            if photon_wvl <= wvl_thresholds[type_index]

                if parent_type == "H2O"
                    tsh_energies = (
                        5.113 , # "H2O + γ -> OH(X2π) + H"
                        7.05, # "H2O + γ -> O + H2"
                        9.12, # "H2O + γ -> OH(A2Σ+) + H"
                        (9.12, 9.54), # "H2O + γ -> O + H + H"
                        12.60, # " H2O + γ -> H2O(+) + e-"
                        (12.60, 18.11), # " H2O + γ -> H + OH(+) + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
                        (12.60, 18.65), # " H2O + γ -> H2 + O(+) + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
                        (12.60, 18.72) # " H2O + γ -> H(+) + OH + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
                    )
                    species_names = (
                        ("H2O", "OH(X2π)", "H"),  # "H2O + γ -> OH(X2π) + H"
                        ("H2O", "O(3P)", "H2"),  # "H2O + γ -> O + H2"
                        ("H2O", "OH(A2Σ+)", "H"), # "H2O + γ -> OH(A2Σ+) + H"
                        (("H2O", "OH(A2Σ+)", "H"), ("OH(A2Σ+)", "O", "H")),          # "H2O + γ -> O + H + H"
                        ("H2O", "H2O(+)", "e(-)"),     # " H2O + γ -> H2O(+) + e-"
                        (("H2O", "H2O(+)", "e(-)"), ("H2O(+)", "OH(+)", "H")),  # " H2O + γ -> H + OH(+) + e-"
                        (("H2O", "H2O(+)", "e(-)"), ("H2O(+)", "O(+)", "H2")),  # " H2O + γ -> H2 + O(+) + e-"
                        (("H2O", "H2O(+)", "e(-)"), ("H2O(+)", "OH", "H(+)"))   # " H2O + γ -> H(+) + OH + e-"
                    )
                    reaction_probabilities = (
                        (0.99, 0.70, 0.70, 0.00, 0.00, 0.00),  # "H2O + γ -> OH(X2π) + H"
                        (0.01, 0.10, 0.10, 0.00, 0.00, 0.00),  # "H2O + γ -> O + H2"
                        (0.00, 0.08, 0.08, 0.00, 0.00, 0.00),  # "H2O + γ -> OH(A2Σ+) + H"
                        (0.00, 0.12, 0.12, 0.00, 0.00, 0.00),  # "H2O + γ -> O + H + H"
                        (0.00, 0.00, 0.00, 1.00, 0.60, 0.63),  # " H2O + γ -> H2O(+) + e-"
                        (0.00, 0.00, 0.00, 0.00, 0.40, 0.31),  # " H2O + γ -> H + OH(+) + e-"
                        (0.00, 0.00, 0.00, 0.00, 0.00, 0.043), # " H2O + γ -> H2 + O(+) + e-"
                        (0.00, 0.00, 0.00, 0.00, 0.00, 0.017)  # " H2O + γ -> H(+) + OH + e-"
                    )
                    reaction_names = (
                        "H2O-PD1", # "H2O + γ -> OH(X2π) + H"
                        "H2O-PD2", # "H2O + γ -> O + H2"
                        "H2O-PD3", # "H2O + γ -> OH(A2Σ+) + H"
                        "H2O-PD4", # "H2O + γ -> O + H + H"
                        "H2O-PI1", # " H2O + γ -> H2O(+) + e-"
                        "H2O-PI2", # " H2O + γ -> H + OH(+) + e-"
                        "H2O-PI3", # " H2O + γ -> H2 + O(+) + e-"
                        "H2O-PI4"  # " H2O + γ -> H(+) + OH + e-"
                    )
                    reaction_types = (
                        "SPD", # "H2O + γ -> OH(X2π) + H"
                        "SPD", # "H2O + γ -> O + H2"
                        "SPD", # "H2O + γ -> OH(A2Σ+) + H"
                        "DPD", # "H2O + γ -> O + H + H"
                        "SPI", # " H2O + γ -> H2O(+) + e-"
                        "DiPI", # " H2O + γ -> H + OH(+) + e-"
                        "DiPI", # " H2O + γ -> H2 + O(+) + e-"
                        "DiPI"  # " H2O + γ -> H(+) + OH + e-"
                    )

                    wavelength_range = ((1357.0, 1860.0), (1208.0, 1220.0), (984.0, 1357.0), (684.4, 984.0), (662.0, 684.4), (0.0, 662.0))

                elseif parent_type == "OH"
                    tsh_energies = (
                        4.39, # " OH(A2Σ+) (v'=3) -> O(3P) + H"
                        4.39, # " OH(A2Σ+) (v'=2) -> O(3P) + H"
                        4.39, # " OH(1Σ+) -> O(3P) + H"
                        2.42, # " OH(12Δ/22Π) -> O(1D) + H"
                        0.20, # " OH(B2Σ) -> O(1S) + H"
                        4.39, # " OH(D2Σ) -> O(3P) + H"
                        13.36, # " OH() -> OH(+) + e-"
                    )
                    species_names = (
                        ("OH(A2Σ+, v'=3)", "O(3P)", "H"), # " OH(A2Σ+) (v'=3) -> O(3P) + H"
                        ("OH(A2Σ+, v'=2)", "O(3P)", "H"), # " OH(A2Σ+) (v'=2) -> O(3P) + H"
                        ("OH(1Σ+)", "O(3P)", "H"), # " OH(1Σ+) -> O(3P) + H"
                        ("OH(12Δ/22Π)", "O(1D)", "H"), # " OH(12Δ/22Π) -> O(1D) + H"
                        ("OH(B2Σ)", "O(1S)", "H"), # " OH(B2Σ) -> O(1S) + H"
                        ("OH(D2Σ)", "O(3P)", "H"), # " OH(D2Σ) -> O(3P) + H"
                        ("OH", "OH(+)", "e(-)"),   # " OH() -> OH(+) + e-"
                    )
                    reaction_probabilities = (
                        (1.00, 0.00, 0.00, 0.00, 0.00, 0.00), # " OH(A2Σ+) (v'=3) -> O(3P) + H"
                        (0.00, 1.00, 0.00, 0.00, 0.00, 0.00), # " OH(A2Σ+) (v'=2) -> O(3P) + H"
                        (0.00, 0.00, 1.00, 0.00, 0.00, 0.00), # " OH(1Σ+) -> O(3P) + H"
                        (0.00, 0.00, 0.00, 0.89, 0.00, 0.00), # " OH(12Δ/22Π) -> O(1D) + H"
                        (0.00, 0.00, 0.00, 0.11, 0.00, 0.00), # " OH(B2Σ) -> O(1S) + H"
                        (0.00, 0.00, 0.00, 0.00, 1.00, 0.00), # " OH(D2Σ) -> O(3P) + H"
                        (0.00, 0.00, 0.00, 0.00, 0.00, 1.00), # " OH() -> OH(+) + e-"
                    )
                    reaction_names = (
                        "OH-PD1", # " OH(A2Σ+) (v'=3) -> O(3P) + H"
                        "OH-PD2", # " OH(A2Σ+) (v'=2) -> O(3P) + H"
                        "OH-PD3", # " OH(1Σ+) -> O(3P) + H"
                        "OH-PD4", # " OH(12Δ/22Π) -> O(1D) + H"
                        "OH-PD5", # " OH(B2Σ) -> O(1S) + H"
                        "OH-PD6", # " OH(D2Σ) -> O(3P) + H"
                        "OH-PI", # " OH() -> OH(+) + e-"
                    )
                    reaction_types = (
                        "SPD", # " OH(A2Σ+) (v'=3) -> O(3P) + H"
                        "SPD", # " OH(A2Σ+) (v'=2) -> O(3P) + H"
                        "SPD", # " OH(1Σ+) -> O(3P) + H"
                        "SPD", # " OH(12Δ/22Π) -> O(1D) + H"
                        "SPD", # " OH(B2Σ) -> O(1S) + H"
                        "SPD", # " OH(D2Σ) -> O(3P) + H"
                        "SPI", # " OH() -> OH(+) + e-"
                    )
                    wavelength_range = ((2439.0,2460.0), (2150.0,2170.0), (1400.0, 1800.0), (1208.0, 1220.0), (1208.0, 1220.0), (928.0, 1000.0), (0.0, 928.0))
                    # wavelength_range = ((2450.0,2460.0), (2140.0,2200.0), (1400.0, 1800.0), (1208.0, 1220.0), (928.0, 1200.0), (0.0, 928.0))
                
                elseif parent_type == "H2"
                    tsh_energies = (
                        4.48, # " H2 -> H(1S) + H(1S)"
                        14.68, # " H2 -> H(2S,2P) + H(2S,2P)"
                        15.43, # " H2 -> H2(+) + e-"
                        (15.43, 17.82) # " H2 -> H(+) + H + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
                    )
                    species_names = (
                        ("H2", "H", "H"), # " H2 -> H(1S) + H(1S)"
                        ("H2", "H", "H"), # " H2 -> H(2S,2P) + H(2S,2P)"
                        ("H2", "H2(+)", "e(-)"),   # " H2 -> H2(+) + e-"
                        (("H2", "H2(+)", "e(-)"), ("H2(+)", "H(+)", "H"))  # " H2 -> H(+) + H + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
                    )
                    reaction_probabilities = (
                        (1.00, 0.00, 0.00, 0.00), # " H2 -> H(1S) + H(1S)"
                        (0.00, 1.00, 0.00, 0.00), # " H2 -> H(2S,2P) + H(2S,2P)"
                        (0.00, 0.00, 1.00, 0.60), # " H2 -> H2(+) + e-"
                        (0.00, 0.00, 0.00, 0.40)  # " H2 -> H(+) + H + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
                    )
                    reaction_names = (
                        "H2-PD1", # " H2 -> H(1S) + H(1S)"
                        "H2-PD2", # " H2 -> H(2S,2P) + H(2S,2P)"
                        "H2-PI1", # " H2 -> H2(+) + e-"
                        "H2-PI2"  # " H2 -> H(+) + H + e-"
                    )
                    reaction_types = (
                        "SPD", # " H2 -> H(1S) + H(1S)"
                        "SPD", # " H2 -> H(2S,2P) + H(2S,2P)"
                        "SPI", # " H2 -> H2(+) + e-"
                        "DiPI"  # " H2 -> H(+) + H + e-"
                    )
                    wavelength_range = ((844.8, 2768.9), (803.7, 844.8), (695.8, 803.7), (0, 695.8))
                
                elseif parent_type == "H"
                    tsh_energies = (13.60) # "H -> H(+) + e-"
                    species_names = (("H", "H(+)", "e(-)"),)
                    reaction_probabilities = ((1.00),)
                    reaction_names = ("H-PI1",)
                    reaction_types = ("SPI",)
                    wavelength_range = ((0, 911.75),)
                
                elseif parent_type == "H(-)"
                    tsh_energies = (0.75,  # "H(-) -> H + e-"
                                    (0.75, 14.63) # "H(-) -> H(+) + 2e-"
                                    ) 
                    species_names = (("H(-)", "H", "e(-)"),
                                     (("H(-)", "H", "e(-)"), ("H", "H(+)", "e(-)"))
                                     )
                    reaction_probabilities = ((1.00, 0.99),
                                              (0.00, 0.01)
                                              )
                    reaction_names = ("H(-)-PI1",
                                      "H(-)-PI2"
                                      )
                    reaction_types = ("SPI",
                                      "DPI" # Double photoionization
                                      ) 
                    wavelength_range = ((847.36, 16439.0), #Lykke et al 1992 and  Huebner 2015
                                        (0, 847.36)
                                        ) 

                
                elseif parent_type == "HO2"
                    tsh_energies = (2.82) # "HO2 -> OH + O"
                    species_names = (("HO2", "OH", "O"),)
                    reaction_probabilities = ((1.00),)
                    reaction_names = ("HO2-PD1",)
                    reaction_types = ("SPD",)
                    wavelength_range = ((0, 4395.0),)
                
                elseif parent_type == "H2O2"
                    tsh_energies = (2.15) # "H -> H(+) + e-"
                    species_names = (("H2O2", "OH", "OH"),)
                    reaction_probabilities = ((1.00),)
                    reaction_names = ("H2O2-PD1",)
                    reaction_types = ("SPD",)
                    wavelength_range = ((0, 5765.0),)

                elseif parent_type == "He"
                    tsh_energies = (24.59) # "He -> He(+) + e-"
                    species_names = (("He", "He(+)", "e(-)"),)
                    reaction_probabilities = ((1.00),)
                    reaction_names = ("He-PI1",)
                    reaction_types = ("SPI",)
                    wavelength_range = ((0, 504.27),)

                elseif parent_type == "Ne"
                    tsh_energies = (21.57) # "Ne -> Ne(+) + e-"
                    species_names = (("Ne", "Ne(+)", "e(-)"),)
                    reaction_probabilities = ((1.00),)
                    reaction_names = ("Ne-PI1",)
                    reaction_types = ("SPI",)
                    wavelength_range = ((0, 574.94),)

                end
            

                for (index, range) in enumerate(wavelength_range)

                    if range[1] < photon_wvl <= range[2]  

                        wvl_range = range
                        
                        # Determine occuring reaction
                        weights = [p[index] for p in reaction_probabilities]
                        reaction_index = rand(Categorical(weights))

                        present_reaction = reaction_types[reaction_index]
                        reaction_name = reaction_names[reaction_index]

                        # Simple photodissociations
                        if present_reaction == "SPD"

                            product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), species_names[reaction_index])

                            if (photon_energy/1.602e-19 >= tsh_energies[reaction_index]) # To avoid problems with H2O-PD2
                                reaction = SimplePhotodissociation.PhotoReaction(tsh_energies[reaction_index], parent_velocity, sun_tuple, nothing, nothing, species_names[reaction_index], false)
                                final_speeds_light, final_speeds_heavy = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)
                                product_types = [product_names[2], product_names[3]]
                                product_velocities = [final_speeds_heavy, final_speeds_light]
                            end

                        # H2O Double photodissociation
                        elseif present_reaction == "DPD"

                            product_names = map(s -> replace(s, r"\((?!\+).*?\)" => ""), species_names[reaction_index])

                            reaction = SimplePhotodissociation.PhotoReaction(tsh_energies[reaction_index][1], parent_velocity, sun_tuple, nothing, nothing, species_names[reaction_index][1], false)
                            final_speeds_Hf, final_speeds_OH = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)

                            reaction = SimplePhotodissociation.PhotoReaction(tsh_energies[reaction_index][2], final_speeds_OH[end], sun_tuple, nothing, nothing, species_names[reaction_index][2], false)
                            final_speeds_Hs, final_speeds_O = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)
                            
                            product_types = ["O", "H", "H"]
                            product_velocities = [final_speeds_O, final_speeds_Hf, final_speeds_Hs]

                        # Double electron ejection of negative atomic hydrogen, equivalente to a double photoionization process
                        elseif present_reaction == "DPI"

                            reaction = SimplePhotoionisation.PhotoReaction(tsh_energies[reaction_index][1], parent_velocity, sun_tuple, nothing, nothing, species_names[reaction_index][1], false)
                            final_speeds_H = SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)

                            reaction = SimplePhotoionisation.PhotoReaction(tsh_energies[reaction_index][2], final_speeds_H[end], sun_tuple, nothing, nothing, species_names[reaction_index][2], false)
                            final_speeds_Hion = SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)
                            
                            product_types = [species_names[reaction_index][2][2]]
                            product_velocities = [final_speeds_Hion]

                        # Simple photoionisations
                        elseif present_reaction == "SPI"

                            reaction = SimplePhotoionisation.PhotoReaction(tsh_energies[reaction_index], parent_velocity, sun_tuple, nothing, nothing, species_names[reaction_index], false)
                            final_speeds_ion = SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)
                            product_types = [species_names[reaction_index][2]]
                            product_velocities = [final_speeds_ion]

                        # Dissociative photoionisations
                        elseif present_reaction == "DiPI"

                            reaction = SimplePhotoionisation.PhotoReaction(tsh_energies[reaction_index][1], parent_velocity, sun_tuple, nothing, nothing, species_names[reaction_index][1], false)
                            final_speeds_ion = SimplePhotoionisation.simulate_photoionisation(reaction, photon_energy)

                            product_names = map(s -> replace(s, r"\(\+\)" => ""), species_names[reaction_index][2])
                            
                            reaction = SimplePhotodissociation.PhotoReaction(tsh_energies[reaction_index][2], final_speeds_ion[end], sun_tuple, nothing, nothing, product_names, false)
                            final_speeds_light, final_speeds_heavy = SimplePhotodissociation.simulate_photodissociation(reaction, photon_energy)
                            
                            product_velocities = [final_speeds_heavy, final_speeds_light]
                            product_types = [species_names[reaction_index][2][2], species_names[reaction_index][2][3]]
                        end

                        break
                    end
                end
            else
                reaction_occurence = false
            end
        end
    end

    return reaction_occurence, reaction_name, product_types, product_velocities, wvl_range

end

end