<table>
<tr>
<td style="vertical-align: top; padding-right: 30px;">

# EXOSPHID

**EXOSpheric PHotoIonisation and PhotoDissociation**

Julia-based scientific code library to simulate photodissociation and photoionisation processes in **extraterrestrial exospheres**, with particular focus on **lunar applications** and **hydrogen-bearing species**.  

Developed by **Carlos Gómez de Olea Ballester** at **TU Munich, Professorship of Lunar and Planetary Exploration**, supervised by **Prof. Dr. Philipp Reiss** and advised by **Alexander Peschel** (PhD candidate).

</td>
<td style="vertical-align: top;">

<img src="https://github.com/user-attachments/assets/b44495b0-c21e-4c1f-830a-51210dd5b003" alt="EXOSPHID Logo" width="500"/>

</td>
</tr>
</table>


---

## ✨ Overview
This repository provides a modular framework to study **photochemical loss pathways** of neutral and ionic species in surface-bounded exospheres of airless bodies.  
Currently, the code supports **10 exospheric species** (`H₂O`, `OH`, `H₂`, `H`, `H⁻`, `HO₂`, `H₂O₂`, `He`, `Ne`, `Ar`) with rates derived from literature and the thesis work.  

Any chemical species may be added by following the methodology described in the [`methods/`](methods/) subfolder.

- **Core focus**: Reproducing solar-driven destruction pathways of lunar water.  
- **Extendable**: Can be adapted for other bodies (Mercury, Europa, asteroids).  
- **Individual particle approach**: Stochastic simulation of photon–particle interactions.

---

## 📖 Reference

A preliminary version of `EXOSPHID` was presented in: 

> *Photolysis of Lunar Water in the Exosphere and on the Surface*,  
> Carlos Gómez de Olea Ballester, TU Munich, 2025.  
> [Available online](https://mediatum.ub.tum.de/node?id=1784196)  

where the package focused on simplified photochemistry for hydrogen-bearing species in the lunar environment ($H_2, OH, H_2O$).

Further work will appear in **upcoming publications (2026)**.

---

## 📚 Documentation
Extensive documentation is available in the [`EXOSPHID WIKI`](https://github.com/charligolea/EXOSPHID/wiki). This includes information regarding:
- Structure and Logic
- General Photochemistry Theory
- Photochemical Database
- Performance & Benchmark
- Validation
- Testing

---

## 📂 Repository Structure
```
├── benchmark/                 
│   └── benchmark.jl            # Benchmarking scripts
│
├── methods/                   
│   ├── add_new_species.ipynb   # Guide for adding new species
│   └── Helium_quiet.csv        # Example input dataset
│
├── src/
│   ├── EXOSPHID.jl             # Core module entry point
│   │
│   ├── database/               # Data and spectrum modules
│   │   ├── photodatabase.jl
│   │   ├── solar_database.jl
│   │   └── solar_spectrum.jl
│   │
│   ├── species/                # Species-specific definitions
│   │   ├── H2O.jl
│   │   ├── OH.jl
│   │   ├── H2.jl
│   │   ├── H.jl
│   │   ├── H(-).jl
│   │   ├── He.jl
│   │   ├── Ne.jl
│   │   ├── Ar.jl
│   │   ├── HO2.jl
│   │   ├── H2O2.jl
│   │   └── general_construct.jl
│   │
│   ├── main/
│   │   └── photodestruction_main.jl  # Main execution script
│   │
│   ├── photoreactions/         # Reaction logic and definitions
│   │   ├── SimplePhotodissociation.jl
│   │   ├── SimplePhotoionisation.jl
│   │   ├── photodestruction_logic.jl
│   │   ├── photodestruction_logic_multiple.jl
│   │   └── photoreactions.jl
│   │
│   ├── speed_distrbutions/     # Example usage
│   │   └── speed_distrbutions.jl
│   │
│   └── test/                   # Testing and validation
│       ├── run_tests.jl
│       ├── photoreactions/
│       │   ├── SimplePhotodissociation.jl
│       │   └── SimplePhotoionisation.jl
│       ├── database/
│       │   ├── photodatabase.jl
│       │   ├── solar_database.jl
│       │   └── solar_spectrum.jl
│       └── validation/
│           └── validation.jl
│
└── README.md                   # Project overview and documentation
```

---
---

## ⚙️ Installation & Usage

Please refer to the [official guide](https://julialang.org/downloads/) to installing **Julia** on your machine.  
The **EXOSPHID.jl** package can be installed using Julia’s package manager **Pkg**:

```julia
using Pkg
Pkg.add(url="https://github.com/charligolea/EXOSPHID.jl.git")
using EXOSPHID
```

If you want to use the latest **development version**, you can install the package as:

```julia
Pkg.add(url="https://github.com/charligolea/EXOSPHID.jl.git", rev="dev")
```

> ⚠️ **Note:** The development version may contain bugs and is not recommended for general use.

---

## 📌 Example Usage: `photodestruction(solar_activity, dt, parent_type, parent_velocity, sun_tuple)`

This function implements the **core simulation engine** for photodestruction of exospheric species within the EXOSPHID framework.  

After correctly installing the `EXOSPHID` package on your Julia machine, you can test the function as follows:

```julia
photodestruction(0, 1e10, "H2O", 550, nothing)
```

- **Inputs**:
  - `solar_activity` : Ratio of solar activity from 0.0 to 1.0, where 0.0 corresponds to Quiet Sun (QS) conditions, and 1.0 is for the Active Sun (AS)
  - `dt` : Simulation time window [s]
  - `parent_type` : `"H2O"`, `"OH"`, `"H2"`, `"H"`, `"H(-)"`, `"HO2"`, `"H2O2"`, `"He"`, `"Ne"`, `"Ar"`
  - `parent_velocity`: Parent velocity [m/s]. Can be given as scalar, or as 3D array (preferably Tuple).
  - `sun_tuple`: Solar vector. May be given as 3D array (preferably Tuple), or as `nothing`.

- **Outputs**:
  - `reaction_occurence::Bool` : True if reaction occurs
  - `present_reaction::String` : Identifier of triggered reaction (`H2O-PD1`, `OH-PI`, etc.)
  - `product_types::Vector{String}` : Reaction products
  - `product_velocities::Vector{NTuple{3, Float32}}` : Velocity vectors of products [m/s]
  - `wvl_range::Tuple{Float32, Float32}` : Photon wavelength range that triggered the reaction

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


### 📚 Selected References
- Huebner, W. F. (2015). *Photoionization and photodissociation rates in solar and blackbody radiation fields.*  
- Huebner, W. F. (1992). *Solar photo rates for planetary atmospheres and atmospheric pollutants.*  
- Combi, M. R. (2005). *Gas Dynamics and Kinetics in the Cometary Coma: Theory and Observations.*  
- Kronebusch, D. & Berkowitz, J. (1976). *PHOTODISSOCIATIVE IONIZATION IN THE 21-41 eV REGION: O2, N2, CO, NO, CO2, H2O, NH3 AND CH4.*  
- Gómez de Olea Ballester, C. (2025). *Photolysis of Lunar Water in the Exosphere and on the Surface*, TU Munich.  

## ⚠️ Notes
- Currently calibrated for 9 species.  
- Rates and pathways extendable via [`src/database/species`](src/database/species/).  
- Parameters and branching ratios subject to update in future publications.

---

## 🧪 Testing
Testing scripts are provided in [`test/`](test/).  
These scripts include physics-based analysis as well as numerical tests.

---

## ✅ Validation
A validation example is provided in [`validation.jl`](validation/validation.jl).  
This script allows comparisons with published results to ensure reproducibility.

---

## 🧩 Performance
Benchmarking examples are provided in [`benchmark.jl`](benchmark/benchmark.jl).  
These scripts include comparisons with published results to ensure reproducibility.

---

## 🚀 Applications
`EXOSPHID` provides a flexible toolset to model photodestruction reactions, which can be integrated within a larger exospheric modeling framework. One such example is:

- [**Extraterrestrial Exosphere and Surface Simulations (ExESS)**](https://github.com/Smolkaa/ExESS.jl)
  Developed by A. Peschel, TU Munich.  
  *(contact: a.peschel@tum.de)*

`ExESS` is a comprehensive numerical model of the lunar exosphere that simulates the distribution and behaviour of multiple volatile species.  

---

## 📜 Citation
If you use this code, please cite (ADD JOSS PAPER HERE WHEN AVAILABLE):

```
Gómez de Olea Ballester, C. & Peschel, A. (2026). EXOSPHID: A Julia Package for Simulating Photodissociation and Photoionisation Processes in Solar System Exospheres. 
```

---

## 👥 Authors & Contact

For questions, please contact:  
- Carlos Gómez de Olea Ballester (primary), carlos.olea@tum.de  
- Alexander Peschel, a.peschel@tum.de
