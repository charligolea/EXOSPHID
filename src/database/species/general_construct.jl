"""
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# GENERAL SPECIES OBJECT
─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
"""

"""
    Species:

A general data structure for defining the photochemical database of a given species.
Each Species instance collects all relevant parameters for modeling photodestruction and photoionisation channels, as well as their probabilities across wavelength intervals.

# INPUTS

- `wvl_threshold::Float32` -> Threshold wavelength (in Å) above which photodestruction is negligible.

- `quiet_rate::Float32` -> Photodestruction rate under quiet Sun conditions (s⁻¹).
- `active_rate::Float32` -> Photodestruction rate under active Sun conditions (s⁻¹).

- `tsh_energies::Tuple` -> Threshold energies (in eV) for each photoprocess channel. The structure varies depending on the reaction type:
    * SPD: single Float32 (dissociation energy).
    * SPI: single Float32 (ionisation threshold).
    * DPD, DPI, DiPI: 2D tuples of Float32 (first-step and full-process energies).

- `species_names::Tuple` -> Parent and product species for each channel. The structure varies depending on the reaction type:
    * SPD/SPI: 3-element tuple (parent, product1, product2).
    * DPD/DPI/DiPI: nested tuples of 3-element tuples describing multi-step processes.
Electronic states may and should be included, e.g. "OH(X²Π)".

- `reaction_probabilities::Tuple` -> Branching ratios (Float32) per wavelength range.
For a given wavelength, probabilities across all channels must sum to 1.

- `reaction_names::Tuple` -> Identifiers for each reaction channel (custom labels for later analysis).

- `reaction_types::Tuple` -> Interaction class, one of:
    * "SPD", "DPD", "SPI", "DPI", "DiPI"
corresponding respectively to single/double photodissociation and ionisation processes.

- `wavelength_range::NTuple{2, Float32}`` -> Photon wavelength intervals (in Å) for which each reaction is active.
Each entry is a 2-tuple of Float32 defining the lower and upper wavelength bounds.

# EXAMPLE:
H2O = Species(
    `wvl_thr_H2O,
    qr_H2O,
    ar_H2O,
    tsh_H2O,
    sn_H2O,
    rp_H2O,
    rn_H2O,
    rt_H2O,
    wv_H2O`
)

# TESTING:
- Tests are available to make sure that the photochemical database is correctly built. See run_tests.jl
"""
struct Species
    wvl_threshold::Float32
    quiet_rate::Float32
    active_rate::Float32
    tsh_energies::Tuple{Vararg{Union{Float32, NTuple{2, Float32}}}}
    species_names:: Tuple{Vararg{Union{NTuple{3, String}, NTuple{2, NTuple{3, String}}}}}
    reaction_probabilities::NTuple{N, Any} where N
    reaction_names::NTuple{N, String} where N
    reaction_types::NTuple{N, String} where N
    wavelength_range::NTuple{N, NTuple{2, Float32}} where N
end

function Species(wvl_threshold::Real, quiet_rate::Real, active_rate::Real,
                 tsh_energies, species_names, reaction_probabilities,
                 reaction_names, reaction_types, wavelength_range)
    Species(Float32(wvl_threshold), Float32(quiet_rate), Float32(active_rate),
            tsh_energies, species_names, reaction_probabilities,
            reaction_names, reaction_types, wavelength_range)
end