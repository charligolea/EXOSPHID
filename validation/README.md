# Validation

This folder contains scripts and tools to validate **EXOSPHID** simulation results, specifically the **exothermic velocities** of dissociation and ionization products.  

---

## Example Usage

The main function for validation is:

```julia
using EXOSPHID

validate_exothermic_velocities(parent_type, solar_activity, num_reactions)
```

### Function Inputs

- `parent_type`  
  Type of parent molecule. Must be one of the **EXOSPHID species**, e.g., `"H2O"`, `"OH"`, `"H2"`.

- `solar_activity`  
  Solar activity level, from `0` (Quiet Sun) to `1` (Active Sun).  For validation against literature values, always use Quiet Sun (`0`).

- `num_reactions`  
  Number of photon interactions to simulate for **every possible reaction**.  

### Specific Example

```julia
# Validate H2O under quiet Sun conditions with 1e6 photon interactions
validate_exothermic_velocities("H2O", 0, 1_000_000)
```

This will simulate the exothermic velocities for `"H2O"` molecules under quiet Sun conditions.

---

### Validation flow

The standard **EXOSPHID** execution relies on the `photodestruction` function. This function first determines the **current reaction** that is occurring and then calls the [`photodestruction_logic.jl`](https://github.com/charligolea/EXOSPHID/tree/main/src/photoreactions/photodestruction_logic.jl) script **once**. Thus, each execution of `photodestruction` simulates a **single photon interaction**.

For the **validation script**, the procedure is slightly different. For a given number of reactions n (typically large, e.g., 100,000 to 1,000,000), we simulate n imaginary photons with energies sampled from the **normalized spectra**. For all of these photons, we iterate through all **possible photoreactions** for the given parent species and call the [`photodestruction_logic_multiple.jl`](https://github.com/charligolea/EXOSPHID/tree/main/src/photoreactions/photodestruction_logic_multiple.jl) script.  

This script recursively simulates each photoreaction **n times**. At the end, it prints **statistics for every specific reaction**, showing the outcomes after simulating it n times. 


### Notes

- Validation is performed by comparing **EXOSPHID estimated velocities** with **literature values**.  
- Literature references are currently available only for `H2O`, `OH`, and `H2`.  
- For more information on the validation methodology, see the [EXOSPHID Wiki: Validation Chapter](https://github.com/charligolea/EXOSPHID/wiki/Validation).

---
