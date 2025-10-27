# Validation

This folder contains scripts and tools to validate **EXOSPHID** simulation results, specifically the **exothermic velocities** of dissociation and ionization products.  

---

## Example Usage

The main function for validation is:

```julia
using EXOSPHID

validate_exothermic_velocities(parent_type, solar_activity, num_reactions)