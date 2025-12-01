# How to Add New Species to `EXOSPHID`

This folder contains scripts and tools to add new exospheric species to the **EXOSPHID** logic.

## 1. In [`photodatabase.jl`](../src/database/photodatabase.jl)

1. Add the new exospheric species to the `exosphid_species` array.  
2. Add necessary parent and daughter masses to `mass_species` and `mass_dict` (in amu).  
3. Include logic for the new species in:
   - `get_photodestruction_rates()`  
   - `get_vibrorotational_energy()`  
   - If needed, `get_electronic_energy_predis()` for electronic states.  

---

## 2. Obtain Wavelength-Dependent Absorption Data

1. Find the wavelength-dependent absorption cross-section for the new species.  
2. Follow the steps in [`add_new_species.ipynb`](add_new_species.ipynb) to generate **standard** and **normalized solar fluxes**.

---

## 3. [`solar_database.jl`](../src/database/solar_database.jl)

1. Add the fluxes obtained in Step 2.  
2. Include the new species logic in:  
   - `get_standard_fluxes()`  
   - `get_normalized_fluxes()`

---


## 4. [`src/database/species`](../src/database/species)

1. Create a new file with the **exact name** of the species as in `exosphid_species`.  
2. Add all necessary arrays for the database, including:
   - Species Names  
   - Reaction Types  
   - Reaction Names  
   - Wavelength Thresholds  
   - Wavelength Ranges  
   - Energy Thresholds  
   - Reaction Probabilities  
   - Photodestruction Rates for Quiet and Active Sun conditions  
3. Reference examples in the same folder or check the [`general_construct.jl`](../src/database/species/general_construct.jl) for the `Struct Species` description:

```julia
using EXOSPHID
EXOSPHID.Species
```

---

## 5. [`validation.jl`](../validation/validation.jl) 

1. Add a new key to the velocities dictionary for the parent species.
2. Add example parent velocities based on thermal velocity calculations (see [WIKI](https://github.com/charligolea/EXOSPHID/wiki/Validation)) or your own method.

---

## 6.  Update the Main README.md

1. Add the new species name wherever species are mentioned.
2. Ensure the new species file is included in the EXOSPHID project tree.