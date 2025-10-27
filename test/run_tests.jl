using Test
if !isdefined(Main, :EXOSPHID)
    include(joinpath(@__DIR__, "..", "src", "EXOSPHID.jl"))
    using .EXOSPHID
end

# Begin EXOSPHID package test

@testset verbose=true "EXOSPHID ............................." begin
    
    @testset verbose=true "PHOTODATABASE ............................." begin
        include(joinpath(@__DIR__, "database", "photodatabase.jl"))
    end

    @testset verbose=true "SOLAR DATABASE ............................." begin
        include(joinpath(@__DIR__, "database", "solar_spectrum.jl"))
        include(joinpath(@__DIR__, "database", "solar_database.jl"))
    end

    @testset verbose=true "PHOTOREACTIONS ............................." begin
        include(joinpath(@__DIR__, "photoreactions", "SimplePhotodissociation.jl"))
        include(joinpath(@__DIR__, "photoreactions", "SimplePhotoionisation.jl"))
    end

end