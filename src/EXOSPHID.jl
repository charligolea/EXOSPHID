module EXOSPHID

# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# LOAD PACKAGES
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

using Random
using LinearAlgebra
using Distributions
using Statistics
using DataFrames


# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# LOAD SOURCES
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "main", "photodestruction_main.jl"))
include(joinpath(@__DIR__, "database", "species", "general_construct.jl"))
include(joinpath(@__DIR__, "database", "photodatabase.jl"))
include(joinpath(@__DIR__, "database", "solar_database.jl"))
include(joinpath(@__DIR__, "database", "solar_spectrum.jl"))

for parent in exosphid_species
    include(joinpath(@__DIR__, "database", "species", "$parent.jl"))
end

include(joinpath(@__DIR__, "photoreactions", "photoreactions.jl"))
include(joinpath(@__DIR__, "photoreactions", "SimplePhotodissociation.jl"))
include(joinpath(@__DIR__, "photoreactions", "SimplePhotoionisation.jl"))
include(joinpath(@__DIR__, "photoreactions", "photodestruction_logic.jl"))
include(joinpath(@__DIR__, "photoreactions", "photodestruction_logic_multiple.jl"))
include(joinpath(@__DIR__, "..", "validation", "validation.jl"))


"""
    EXOSPHID.jl

## EXOSpheric PhotoIonisation & PhotoDissociation

Scientific code library for simulating photochemical loss pathways in exosphere environments of airless
bodies in the solar system. Note that the package is under development and subject to
frequent changes. For questions, please contact the C. Gómez de Olea B. via email.

The documentation and manual can be found in the [GitHub wiki](https://github.com/charligolea/EXOSPHID/wiki).

_Author: C. Gómez de Olea Ballester (carloas.olea@tum.de)_

## Installation & Usage

Please refer to the [official guide](https://julialang.org/downloads/platform/) to
installing Julia on your machine. The `EXOSPHID.jl` package can be installed using Julia's
package-manager `Pkg`:
```julia
using Pkg
Pkg.add("https://github.com/charligolea/EXOSPHID.jl.git")
using EXOSPHID
```
If you want to use the latest development version, you can install the package as:
```julia
Pkg.add("https://github.com/charligolea/EXOSPHID.jl.git#dev")
```
"""
EXOSPHID
end
