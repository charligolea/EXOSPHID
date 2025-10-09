function get_photodestruction_rates(species::String, solar_activity::Float32)
    """
    Return the photodestruction rate `k` [1/s] 
    for a given `species` and `solar_activity` level (0 = quiet, 1 = active).
    """

    @assert species in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne") "Invalid parent species: $parent_type"
    @assert 0.0 <= solar_activity <= 1.0 "Solar activity must be in (0,1)!"
    
    if species == "H2O"
        quiet_rate = 12.94f-6   # Gomez de Olea (2026)
        active_rate = 22.25f-6
    elseif species == "OH"
        quiet_rate = 22.45f-6   # Gomez de Olea (2025)
        active_rate = 32.56f-6
    elseif species == "H2"
        quiet_rate = 14.24f-6   # Gomez de Olea (2025)
        active_rate = 32.83f-6
    elseif species == "H"
        quiet_rate = 7.00f-8    # Gomez de Olea (2026)
        active_rate = 1.68f-7
    elseif species == "H(-)"
        quiet_rate = 14.24f0      # Gomez de Olea (2026)
        active_rate = 14.24f0
    elseif species == "HO2"
        quiet_rate = 6.59f-3    # Gomez de Olea (2026)
        active_rate = 6.59f-3
    elseif species == "H2O2"
        quiet_rate = 1.43f-4    # Gomez de Olea (2026)
        active_rate = 1.43f-4
    elseif species == "He"
        quiet_rate = 5.57f-8    # Gomez de Olea (2026)
        active_rate = 16.80f-8
    elseif species == "Ne"
        quiet_rate = 1.72f-7    # Gomez de Olea (2026)
        active_rate = 5.74f-7
    end

    k = (1 - solar_activity) * quiet_rate + solar_activity * active_rate

    return k
end

get_photodestruction_rates(species::String, solar_activity::Real) = get_photodestruction_rates(species, Float32(solar_activity))

function is_photoreaction_occuring(k::Float32, dt::Float32)
    @assert 0.0 <= k "Photodestruction rate must be positive!"
    @assert 0.0 <= dt "dt must be positive!"
    p = 1-exp(-k*dt)
    return rand() < p
end

is_photoreaction_occuring(k::Real, dt::Real) = is_photoreaction_occuring(Float32(k), Float32(dt))

function get_wvl_threshold(parent_type::String)
    @assert parent_type in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne") "Invalid parent species: $parent_type"
    # In Armstrong
    if parent_type == "H2O"
        wvl_threshold =  1860.0f0
    elseif parent_type == "OH"
        wvl_threshold = 2460.0f0
    elseif parent_type == "H2"
        wvl_threshold = 2768.0f0
    elseif parent_type == "H"
        wvl_threshold = 911.75f0
    elseif parent_type == "H(-)"
        wvl_threshold = 16439.0f0
    elseif parent_type == "HO2"
        wvl_threshold = 4395.0f0
    elseif parent_type == "H2O2"
        wvl_threshold = 5765.0f0
    elseif parent_type == "He"
        wvl_threshold = 504.27f0
    elseif parent_type == "Ne"
        wvl_threshold = 574.94f0
    end

    @assert wvl_threshold>0.0 "Wavelength threshold cannot be negative"

    return wvl_threshold

end


function get_species_photochemical_info(parent_type::String)

    @assert parent_type in ("H2O", "OH", "H2", "H", "H(-)", "HO2", "H2O2", "He", "Ne") "Invalid parent species: $parent_type"

    if parent_type == "H2O"
        tsh_energies = (
            5.113f0 , # "H2O + γ -> OH(X2π) + H"
            6.98f0, # "H2O + γ -> O + H2" # New value from Sumin Yan et al. 2021
            9.12f0, # "H2O + γ -> OH(A2Σ+) + H"
            (9.12f0, 9.54f0 + 0.25f0), # "H2O + γ -> O + H + H"
            12.60f0, # " H2O + γ -> H2O(+) + e-"
            (12.60f0, 18.11f0), # " H2O + γ -> H + OH(+) + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
            (12.60f0, 18.65f0), # " H2O + γ -> H2 + O(+) + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
            (12.60f0, 18.72f0) # " H2O + γ -> H(+) + OH + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
        )
        species_names = (
            ("H2O", "OH(X2π)", "H"),  # "H2O + γ -> OH(X2π) + H"
            ("H2O", "O(3P)", "H2"),  # "H2O + γ -> O + H2"
            ("H2O", "OH(A2Σ+)", "H"), # "H2O + γ -> OH(A2Σ+) + H"
            (("H2O", "OH(A2Σ+)", "H"), ("OH(DPD)", "O", "H")),          # "H2O + γ -> O + H + H"
            ("H2O", "H2O(+)", "e(-)"),     # " H2O + γ -> H2O(+) + e-"
            (("H2O", "H2O(+)", "e(-)"), ("H2O(+)", "OH(+)", "H")),  # " H2O + γ -> H + OH(+) + e-"
            (("H2O", "H2O(+)", "e(-)"), ("H2O(+)", "O(+)", "H2")),  # " H2O + γ -> H2 + O(+) + e-"
            (("H2O", "H2O(+)", "e(-)"), ("H2O(+)", "OH", "H(+)"))   # " H2O + γ -> H(+) + OH + e-"
        )
        reaction_probabilities = (
            (1.00f0, 0.99f0, 0.70f0, 0.70f0, 0.00f0, 0.00f0, 0.00f0),  # "H2O + γ -> OH(X2π) + H"
            (0.00f0, 0.01f0, 0.10f0, 0.10f0, 0.00f0, 0.00f0, 0.00f0),  # "H2O + γ -> O + H2"
            (0.00f0, 0.00f0, 0.08f0, 0.08f0, 0.00f0, 0.00f0, 0.00f0),  # "H2O + γ -> OH(A2Σ+) + H"
            (0.00f0, 0.00f0, 0.12f0, 0.12f0, 0.00f0, 0.00f0, 0.00f0),  # "H2O + γ -> O + H + H"
            (0.00f0, 0.00f0, 0.00f0, 0.00f0, 1.00f0, 0.60f0, 0.63f0),  # " H2O + γ -> H2O(+) + e-"
            (0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.40f0, 0.31f0),  # " H2O + γ -> H + OH(+) + e-"
            (0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.043f0), # " H2O + γ -> H2 + O(+) + e-"
            (0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.017f0)  # " H2O + γ -> H(+) + OH + e-"
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

        wavelength_range = ((1760.0f0, 1860.0f0), (1357.0f0, 1760.0f0), (1208.0f0, 1220.0f0), (984.0f0, 1357.0f0), (684.4f0, 984.0f0), (662.0f0, 684.4f0), (0.0f0, 662.0f0))

    elseif parent_type == "OH"
        tsh_energies = (
            4.39f0, # " OH(A2Σ+) (v'=3) -> O(3P) + H"
            4.39f0, # " OH(A2Σ+) (v'=2) -> O(3P) + H"
            4.39f0, # " OH(1Σ+) -> O(3P) + H"
            2.42f0, # " OH(12Δ/22Π) -> O(1D) + H"
            0.20f0, # " OH(B2Σ) -> O(1S) + H"
            4.39f0, # " OH(D2Σ) -> O(3P) + H"
            13.36f0, # " OH() -> OH(+) + e-"
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
            (1.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0), # " OH(A2Σ+) (v'=3) -> O(3P) + H"
            (0.00f0, 1.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0), # " OH(A2Σ+) (v'=2) -> O(3P) + H"
            (0.00f0, 0.00f0, 1.00f0, 0.00f0, 0.00f0, 0.00f0), # " OH(1Σ+) -> O(3P) + H"
            (0.00f0, 0.00f0, 0.00f0, 0.89f0, 0.00f0, 0.00f0), # " OH(12Δ/22Π) -> O(1D) + H"
            (0.00f0, 0.00f0, 0.00f0, 0.11f0, 0.00f0, 0.00f0), # " OH(B2Σ) -> O(1S) + H"
            (0.00f0, 0.00f0, 0.00f0, 0.00f0, 1.00f0, 0.00f0), # " OH(D2Σ) -> O(3P) + H"
            (0.00f0, 0.00f0, 0.00f0, 0.00f0, 0.00f0, 1.00f0), # " OH() -> OH(+) + e-"
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
        wavelength_range = ((2439.0f0,2460.0f0), (2150.0f0,2170.0f0), (1400.0f0, 1800.0f0), (1208.0f0, 1220.0f0), (928.0f0, 1000.0f0), (0.0f0, 928.0f0))
        # wavelength_range = ((2450.0,2460.0), (2140.0,2200.0), (1400.0, 1800.0), (1208.0, 1220.0), (928.0, 1200.0), (0.0, 928.0))
    
    elseif parent_type == "H2"
        tsh_energies = (
            4.48f0, # " H2 -> H(1S) + H(1S)"
            14.68f0, # " H2 -> H(2S,2P) + H(2S,2P)"
            15.43f0, # " H2 -> H2(+) + e-"
            (15.43f0, 17.82f0) # " H2 -> H(+) + H + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
        )
        species_names = (
            ("H2", "H", "H"), # " H2 -> H(1S) + H(1S)"
            ("H2", "H", "H"), # " H2 -> H(2S,2P) + H(2S,2P)"
            ("H2", "H2(+)", "e(-)"),   # " H2 -> H2(+) + e-"
            (("H2", "H2(+)", "e(-)"), ("H2(+)", "H(+)", "H"))  # " H2 -> H(+) + H + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
        )
        reaction_probabilities = (
            (1.00f0, 0.00f0, 0.00f0, 0.00f0), # " H2 -> H(1S) + H(1S)"
            (0.00f0, 1.00f0, 0.00f0, 0.00f0), # " H2 -> H(2S,2P) + H(2S,2P)"
            (0.00f0, 0.00f0, 1.00f0, 0.60f0), # " H2 -> H2(+) + e-"
            (0.00f0, 0.00f0, 0.00f0, 0.40f0)  # " H2 -> H(+) + H + e-" (The first is the ionisation energy, the second is the total dissociative ionisation energy)
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
        wavelength_range = ((844.8f0, 2768.9f0), (803.7f0, 844.8f0), (695.8f0, 803.7f0), (0f0, 695.8f0))
    
    elseif parent_type == "H"
        tsh_energies = (13.60f0) # "H -> H(+) + e-"
        species_names = (("H", "H(+)", "e(-)"),)
        reaction_probabilities = ((1.00f0),)
        reaction_names = ("H-PI1",)
        reaction_types = ("SPI",)
        wavelength_range = ((0f0, 911.75f0),)
    
    elseif parent_type == "H(-)"
        tsh_energies = (0.75f0,  # "H(-) -> H + e-"
                        (0.75f0, 14.63f0) # "H(-) -> H(+) + 2e-"
                        ) 
        species_names = (("H(-)", "H", "e(-)"),
                            (("H(-)", "H", "e(-)"), ("H", "H(+)", "e(-)"))
                            )
        reaction_probabilities = ((1.00f0, 0.99f0),
                                    (0.00f0, 0.01f0)
                                    )
        reaction_names = ("H(-)-PI1",
                            "H(-)-PI2"
                            )
        reaction_types = ("SPI",
                            "DPI" # Double photoionization
                            ) 
        wavelength_range = ((847.36f0, 16439.0f0), #Lykke et al 1992 and  Huebner 2015
                            (0f0, 847.36f0)
                            ) 

    
    elseif parent_type == "HO2"
        tsh_energies = (2.82f0) # "HO2 -> OH + O"
        species_names = (("HO2", "OH", "O"),)
        reaction_probabilities = ((1.00f0),)
        reaction_names = ("HO2-PD1",)
        reaction_types = ("SPD",)
        wavelength_range = ((0f0, 4395.0f0),)
    
    elseif parent_type == "H2O2"
        tsh_energies = (2.15f0) # "H -> H(+) + e-"
        species_names = (("H2O2", "OH", "OH"),)
        reaction_probabilities = ((1.00f0),)
        reaction_names = ("H2O2-PD1",)
        reaction_types = ("SPD",)
        wavelength_range = ((0f0, 5765.0f0),)

    elseif parent_type == "He"
        tsh_energies = (24.59f0) # "He -> He(+) + e-"
        species_names = (("He", "He(+)", "e(-)"),)
        reaction_probabilities = ((1.00f0),)
        reaction_names = ("He-PI1",)
        reaction_types = ("SPI",)
        wavelength_range = ((0f0, 504.27f0),)

    elseif parent_type == "Ne"
        tsh_energies = (21.57f0) # "Ne -> Ne(+) + e-"
        species_names = (("Ne", "Ne(+)", "e(-)"),)
        reaction_probabilities = ((1.00f0),)
        reaction_names = ("Ne-PI1",)
        reaction_types = ("SPI",)
        wavelength_range = ((0f0, 574.94f0),)

    end

    return tsh_energies, species_names, reaction_probabilities, reaction_names, reaction_types, wavelength_range
end

get_reaction_index(photo_info, index) = rand(Categorical([p[index] for p in photo_info[3]]))
get_reaction_identifier(photo_info, reaction_index) = photo_info[4][reaction_index]
get_reaction_type(photo_info, reaction_index) = photo_info[5][reaction_index]
get_wavelength_range(photo_info) = photo_info[6]


function get_photoreaction_characteristics(photon_wvl::Real, photochemical_info::NTuple{6, Tuple})

    wavelength_range = get_wavelength_range(photochemical_info)
    present_reaction = nothing
    reaction_name = nothing
    reaction_index = nothing
    wvl_range = nothing

    for (index, range) in enumerate(wavelength_range)
        if range[1] < photon_wvl <= range[2]  
            wvl_range = range

            # Determine occuring reaction
            reaction_index = get_reaction_index(photochemical_info, index)
            present_reaction = get_reaction_type(photochemical_info, reaction_index)
            reaction_name = get_reaction_identifier(photochemical_info, reaction_index)
            break
        end
    end

    if present_reaction !== nothing
        return present_reaction, reaction_name, reaction_index, wvl_range
    else
        return nothing, nothing, nothing, nothing
    end

end