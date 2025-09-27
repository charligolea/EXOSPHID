# EXOSPHID
**EXOSpheric PHotoIonisation and PhotoDissociation**

Julia-based scientific code library to simulate photodissociation and photoionisation processes in **extraterrestrial exospheres**, with particular focus on **lunar applications** and **hydrogen-bearing species**.  
Developed as part of the Master's Thesis of **Carlos Gómez de Olea Ballester** at **TU Munich, Professorship of Lunar and Planetary Exploration**, supervised by **Prof. Dr. Philipp Reiss** and advised by **Alexander Peschel** (PhD candidate).

---

## ✨ Overview
This repository provides a modular framework to study **photochemical loss pathways** of neutral and ionic species in surface-bounded exospheres of airless bodies.  
Currently, the code supports **9 atomic and molecular species** (`H₂O`, `OH`, `H₂`, `H`, `H⁻`, `HO₂`, `H₂O₂`, `He`, `Ne`) with rates derived from literature and the thesis work.  

Any chemical species may be added by following the methodology described in the [`methods/`](methods/) subfolder.

- **Core focus**: Reproducing solar-driven destruction pathways of lunar water.  
- **Extendable**: Can be adapted for other bodies (Mercury, Europa, asteroids).  
- **Monte Carlo approach**: Stochastic simulation of photon–particle interactions.

---

## 📖 Reference
The theory behind this model and its lunar-specific mechanisms is described in:  

> *Photolysis of Lunar Water in the Exosphere and on the Surface*,  
> Carlos Gómez de Olea Ballester, TU Munich, 2025.  
> [Available online](https://mediatum.ub.tum.de/node?id=1784196)  

Further work will appear in **upcoming publications (2026)**.

---

## 📂 Repository Structure
```
├── src/                     # Source code modules
│   ├── photoreactions.jl    # Reaction definitions
│   ├── solar_spectrum.jl    # Solar flux & spectrum functions
│   └── ...
├── methods/                 # Methodological documentation & species expansion guide
├── testing/                 # Testing & benchmarking scripts
│   └── validation.jl
├── photodestruction_main.jl # Main executable script
└── README.md                # This file
```

---

## 🚀 Usage
This repository is designed to be integrated with an exospheric framework such as:

- **Extraterrestrial Exosphere and Surface Simulations (EESS)**  
  Developed by A. Peschel, TU Munich.  
  *(contact: a.peschel@tum.de)*

### 🔧 Main Function: `photodestruction(solar_activity, dt, parent_type, parent_velocity, sun_tuple)`

This script implements the **core simulation engine** for photodestruction of exospheric species within the EXOSPHID framework.  

It provides a full database of photodestruction pathways for the currently supported **9 atomic and molecular species**, including:  
- Threshold energies  
- Branching ratios  
- Wavelength ranges of all known photodestruction mechanisms  

These values are sourced from a range of literature (**Huebner 2015**, Combi 2005, Kronebusch & Berkowitz 1976, Lykke 1992) and expanded by Carlos Gómez de Olea Ballester (2025,2026 TU Munich).  

Photodestruction rates are species-dependent and calculated for both quiet and active Sun conditions. Rates are combined with photon fluxes to determine whether a given parent molecule undergoes photon-induced photodestruction. If a photon interaction occurs, photon energy and wavelength are sampled from **species-dependent solar flux data** implemented in `src/solar_spectrum.jl`.

- **Inputs**:
  - `solar_activity::Float64` : Ratio of solar activity from 0.0 to 1.0, where 0.0 corresponds to Quiet Sun conditions, and 1.0 is for the Active Sun
  - `dt::Float64/Int64` : Simulation time window [s]
  - `parent_type::String` : `"H2O"`, `"OH"`, `"H2"`, `"H"`, `"H(-)"`, `"HO2"`, `"H2O2"`, `"He"`, `"Ne"`
  - `parent_velocity::Tuple{Float64, Float64, Float64}` or `Float64` : Parent velocity [m/s]
  - `sun_tuple::Tuple{Float64, Float64, Float64}` or `Nothing` : Sun–particle geometry

- **Outputs**:
  - `reaction_occurence::Bool` : True if reaction occurs
  - `present_reaction::String` : Identifier of triggered reaction (`H2O-PD1`, `OH-PI`, etc.)
  - `product_types::Array{String}` : Reaction products
  - `product_velocities::Array{Tuple}` : Velocity vectors of products [m/s]
  - `wvl_range::Tuple` : Photon wavelength range that triggered the reaction

---

### 🧪 Supported Photochemical Reaction Types

The model accounts for a range of fundamental processes relevant in exospheric environments:

- **Simple Photodissociation (SPD)**  
  Parent molecule breaks into two fragments after photon absorption.  
  *Example:* `OH + γ → O + H`

- **Double Photodissociation (DPD)**  
  Parent molecule undergoes sequential breaking into three fragments.  
  *Example:* `H2O + γ → O + H + H`

- **Simple Photoionisation (SPI)**  
  Single ionisation with ejection of one electron.  
  *Example:* `H2 + γ → H2(+) + e(-)`

- **Dissociative Photoionisation (DiPI)**  
  Photon absorption both ionises and dissociates the parent.  
  *Example:* `H2O + γ → H + OH(+) + e(-)`

- **Double Photoionisation (DPI)**  
  Double electron ejection (special case for negative hydrogen).  
  *Example:* `H(-) + γ → H(+) + 2e(-)`


### 📚 References
- Huebner, W. F. (2015). *Photoionization and photodissociation rates in solar and blackbody radiation fields.*  
- Huebner, W. F. (1992). *Solar photo rates for planetary atmospheres and atmospheric pollutants.*  
- Combi, M. R. (2005). *Gas Dynamics and Kinetics in the Cometary Coma: Theory and Observations.*  
- Kronebusch, D. & Berkowitz, J. (1976). *PHOTODISSOCIATIVE IONIZATION IN THE 21-41 eV REGION: O2, N2, CO, NO, CO2, H2O, NH3 AND CH4.*  
- Gómez de Olea Ballester, C. (2025). *Photolysis of Lunar Water in the Exosphere and on the Surface*, TU Munich.  

---

## 🔗 Integration
This script is intended for integration with **Extraterrestrial Exosphere and Surface Simulations (EESS)** by A. Peschel.  
It provides reaction rates and Monte Carlo outcomes for use in higher-level **exospheric transport and dynamics models**.

---

## 📌 Example Usage
```julia
include("photodestruction_main.jl")

# Example: Calculate reaction probability for H2O under active Sun
reaction_occurence, reaction_name, products, velocities, wvl_range =
    photodestruction(1.0, 100.0, "H2O", (0.0, 0.0, 0.0), nothing)

println(reaction_name, products, velocities)
```

---

## ⚠️ Notes
- Currently calibrated for 9 species.  
- Rates and pathways extendable via [`src/photoreactions.jl`](../src/photoreactions.jl).  
- Parameters and branching ratios subject to update in future publications.

---

## 🧪 Testing
Validation and benchmarking examples are provided in [`testing/`](testing/).  
These scripts include comparisons with published results to ensure reproducibility.

---

## 📜 Citation
If you use this code, please cite:

```
Gómez de Olea Ballester, C. (2025). Photolysis of Lunar Water in the Exosphere and on the Surface. 
Technical University of Munich, Professorship of Lunar and Planetary Exploration.
Available at: https://mediatum.ub.tum.de/node?id=1784196
```

---

## 👥 Authors & Contact
- **Carlos Gómez de Olea Ballester**, 2025  
  Master’s Thesis, TU Munich  

**Supervised by**: Prof. Dr. Philipp Reiss  
**Advised by**: Alexander Peschel (PhD candidate, TU Munich)

For questions, please contact:  
- Carlos Gómez de Olea Ballester (primary), carlos.olea@tum.de  
- Alexander Peschel, a.peschel@tum.de
