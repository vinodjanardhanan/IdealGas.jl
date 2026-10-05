using IdealGas
using Test

@testset "IdealGas.jl" begin
    
    @testset "Thermoall test " begin        
        thObj = create_thermo(["CH4", "CO", "CO2", "H2", "H2O", "O2"], "lib/therm.dat" )
        cpch4 = IdealGas.cp("CH4", 298.15, thObj)
        hch4 = IdealGas.H("CH4", 298.15, thObj)
        sch4 = IdealGas.S("CH4", 298.15, thObj)        
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
        eh2 = nernst_potential(H2Oxidation(), E0, 1073.15, aH2=aH2, aO2=aO2, aH2O=aH2O)                
        @test 0.9 < eh2 < 1.1        
        eh2 = nernst_h2(E0, 1073.15, pH2=1e5, pO2=0.21e5, pH2O=1.0, water=:liquid)                
        @test 0.9 < eh2 < 1.1        
        eh2 = nernst_h2(E0, 1073.15, pH2=1e5, pO2=0.21e5, pH2O=0.79e5, water=:gas)                
        @test 0.9 < eh2 < 1.1        
        eh2 = nernst_h2(E0, 1073.15, pH2=1e5, pO2=0.21e5, pH2O=0.79e5, water=:vapor)                
        @test 0.9 < eh2 < 1.1        
    end

    @testset "Nernst potential for CO" begin        
        thObj = create_thermo(["CO", "CO2", "O2"], "lib/therm.dat" )
        E0 = E0_CO(thObj,1073.15)       
        aCO = 0.3
        aCO2 = 0.6
        aO2 = 0.21
        eco = nernst_potential(COOxidation(), E0, 1073.15, aCO=aCO, aO2=aO2, aCO2=aCO2)        
        @test 0.9 < eco < 1.1                
        p_std = 1e5
        eco = nernst_co(E0, 1073.15, pCO=1e5, pO2=0.21*p_std, pCO2=0.79*p_std)        
        @test 0.9 < eco < 1.1        
    end

    
end


