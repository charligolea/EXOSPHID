# How to Add New Species to `EXOSPHID`

This folder contains scripts and tools to evaluate the performance of **EXOSPHID**.

---

## Workflow

The current `EXOSPHID` database is mainly suited for studies of the lunar exosphere.  
If you want to **add new species** to the database, please follow the steps in the Jupyter Notebook [`add_new_species.ipynb`](add_new_species.ipynb) and make sure you comply with the following checklist:

---

### 1. In [`photodatabase.jl`](./photodatabase.jl)

1. Add the new atomic/molecular species to the `exosphid_species` array.  
2. Add necessary parent and daughter masses to `mass_species` and `mass_dict`.  
3. Add logic for the new species in `get_species_photochemical_info()`.  
4. Add vibro-rotational energies for all relevant parent and daughter electronic states in `get_vibrorotational_energy()`.  
5. If applicable, do the same for electronic energies in `get_electronic_energy_predis()`.

---

### 2. In [`src/database/species`](./src/database/species)

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

Refer to the [Photochemical Database](https://github.com/charligolea/EXOSPHID/wiki/Photochemical-Database) Wiki chapter for further information.

---

### 3. In [`solar_database.jl`](./src/database/solar_database.jl)

1. Add **standard fluxes** for Quiet and Active Sun, following the existing format.  
2. Repeat the same for **normalized fluxes**.