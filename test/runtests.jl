using IdealGas
using Test

@testset "IdealGas.jl" begin
    retcode = create_thermo(["CH4"], "lib/therm.dat" )
    @test retcode.thermo_all[1].name == "CH4"
end
