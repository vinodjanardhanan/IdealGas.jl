using IdealGas
using Test

@testset "IdealGas.jl" begin
    
    @testset "Thermoall test " begin        
        retcode = create_thermo(["CH4"], "lib/therm.dat" )
        @test retcode.thermo_all[1].name == "CH4"                    
    end

    @testset "STD potential H2" begin
        thObj = create_thermo(["H2", "H2O", "O2"], "lib/therm.dat" )
        E0 = E0_H2(thObj,1073.15)
        @test E0 < 0.977 && E0 > 0.976
    end

    @testset "STD potential CO" begin
        thObj = create_thermo(["CO", "CO2", "O2"], "lib/therm.dat" )
        E0 = E0_CO(thObj,1073.15)
        @test E0 < 0.999 && E0 > 0.98
    end

    @testset "Nernst potential H2-O2 below 100 C" begin
        thObj = create_thermo(["H2", "H2O", "O2"], "lib/therm.dat" )
        E0 = E0_H2(thObj,298.15)
        E_rev = nernst(E0,298.15,pH2=1e5,pO2=0.21e5,pH2O=1e5)                
        @test E_rev > 0 && E_rev < 1.5
    end

    @testset "Nernst potential H2-O2 above 100 C" begin
        thObj = create_thermo(["H2", "H2O", "O2"], "lib/therm.dat" )
        E0 = E0_H2(thObj,1073.15)
        E_rev = nernst(E0,1073.15,pH2=0.5e5,pO2=0.21e5,pH2O=0.5e5)                
        @test E_rev > 0 && E_rev < 1.5
    end

    @testset "Nernst potential CO-CO2 " begin
        thObj = create_thermo(["CO", "CO2", "O2"], "lib/therm.dat" )
        E0 = E0_CO(thObj,1073.15)
        E_rev = nernst_co(E0,1073.15,pCO=0.5e5, pO2=0.21e5,pCO2=0.5e5)                
        @test E_rev > 0 && E_rev < 1.5
    end
end
