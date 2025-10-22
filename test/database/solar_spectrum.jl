include(joinpath(@__DIR__, "src/database/solar_database.jl"))
import solar_database as sd

include(joinpath(@__DIR__, "src/database/solar_spectrum.jl"))
import solar_spectrum as ss

@testset verbose=true "solar_spectrum.jl" begin

end