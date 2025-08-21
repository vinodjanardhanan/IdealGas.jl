using IdealGas
using Test

@testset "IdealGas.jl" begin
    
    @testset "Thermoall test " begin        
        thObj = create_thermo(["CH4", "CO", "CO2", "H2", "H2O", "O2"], "lib/therm.dat" )
        cpch4 = IdealGas.cp("CH4", 298.15, thObj)
        hch4 = IdealGas.H("CH4", 298.15, thObj)
        sch4 = IdealGas.S("CH4", 298.15, thObj)
        println("cp CH4: ", cpch4)
        println("H CH4: ", hch4)
        println("S CH4: ", sch4)
        hall = H_all(thObj, 298.15)
        sall = S_all(thObj, 298.15)
        cpall = cp_all(thObj, 298.15)
        @test hall[1] == hch4
        @test sall[1] == sch4
        @test cpall[1] == cpch4        
    end

    @testset "Nernst potential for H2" begin        
        thObj = create_thermo(["H2", "H2O", "O2"], "lib/therm.dat" )
        E0 = E0_H2(thObj,1073.15)       
        aH2 = 0.3
        aH2O = 0.6
        aO2 = 0.21
        npH2 = NernstH2(E0, 1073.15, aH2, aO2, aH2O)
        eh2 = nerst_potential(npH2)
        @test 0.9 < eh2 < 1.1        
    end

    @testset "Nernst potential for CO" begin        
        thObj = create_thermo(["CO", "CO2", "O2"], "lib/therm.dat" )
        E0 = E0_CO(thObj,1073.15)       
        aCO = 0.3
        aCO2 = 0.6
        aO2 = 0.21
        npCO = NernstCO(E0, 1073.15, aCO, aO2, aCO2)
        eco= nerst_potential(npCO)
        @test 0.9 < eco < 1.1        
    end

end
