"""
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# GENERAL SPECIES OBJECT
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
"""

struct Species
    wvl_threshold::Float32
    quiet_rate::Float32
    active_rate::Float32
    tsh_energies::Tuple
    species_names::Tuple
    reaction_probabilities::Tuple
    reaction_names::Tuple
    reaction_types::Tuple
    wavelength_range::Tuple
end