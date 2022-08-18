module IdealGas
include("Constants.jl")

export create_thermo, cp, H, S, cp_all, cpmix, H_all, Hmix, S_all, Smix, Gmix


abstract type ComponentDefinition end
abstract type ThermoData end

struct Gasphase <: ComponentDefinition    
    species::Array{String,1}
end

struct NASAThermo{T1 <: Float64, T2 <: Integer} <: ThermoData
    name::String
    phase::String
    composition::Dict{String,T2}
    molWt::T1
    ltl::T1
    htl::T1
    cmt::T1
    ltp::Array{T1,1}
    htp::Array{T1,1}
end

#This array will be order as per ig.species
struct SpeciesThermoObj
    molwt::Array{Float64,1}
    thermo_all::Array{NASAThermo}
end
export SpeciesThermoObj


"""
    Tuple of element weights. This is used for the calculation of molecular weights
"""
elementWeight = (    H = 1.00794e-3,  He = 4.002602e-3,  Li = 6.941e-3,  Be = 9.012182e-3,  B = 1.811e-3,
                C = 12.011e-3,  N = 14.00674e-3,  O = 15.9994e-3,  F = 18.9984032e-3,  Ne = 20.1797e-3,
                Na = 22.98977e-3, Mg = 24.3050e-3, Al = 26.98154e-3, Si = 28.0855e-3, P = 30.97376e-3,
                S = 32.066e-3, Cl = 35.4527e-3, Ar = 39.948e-3, K = 39.0983e-3, Ca = 40.078e-3,
                Sc = 44.95591e-3, Ti = 47.88e-3, V = 50.9415e-3, Cr = 51.9961e-3, Mn = 54.9381e-3,
                Fe = 55.847e-3, Co = 58.9332e-3, Ni = 58.69e-3, Cu = 63.546e-3, Zn = 65.39e-3,
                Ga = 69.723e-3, Ge = 72.61e-3, As = 74.92159e-3, Se = 78.96e-3, Br = 79.904e-3,
                Kr = 83.80e-3,  Rb = 85.4678e-3, Sr = 87.62e-3, Y = 88.90585e-3, Zr = 91.224e-3,
                Nb = 92.90638e-3, Mo = 95.94e-3, Tc = 97.9072e-3, Ru = 101.07e-3, Rh = 102.9055e-3,
                Pd = 106.42e-3, Ag = 107.8682e-3,   Cd = 112.411e-3,    In = 114.82e-3, Sn = 118.710e-3,
                Sb = 121.75e-3, Te = 127.6e-3,  I = 126.90447e-3,   Xe = 131.29e-3,   Cs = 132.90543e-3,
                Ba = 137.327e-3, La = 138.9055e-3, Ce = 140.115e-3, Pr = 140.90765e-3, Nd = 144.24e-3,
                Pm = 144.9127e-3, Sm = 150.36e-3,   Eu = 151.965e-3,    Gd = 157.25e-3, Tb = 158.92534e-3,
                Dy = 162.50e-3, Ho = 164.93032e-3,  Er = 167.26e-3, Tm = 168.93421e-3,  Yb = 173.04e-3,
                Lu = 174.967e-3,    Hf = 178.49e-3, Ta = 180.9479e-3,   W = 183.85e-3, Re = 186.207e-3,
                Pt = 195.08e-3, Au = 196.96654e-3,  Hg = 200.59e-3, Tl = 204.3833e-3, Pb = 207.2e-3,
                Bi = 208.98037e-3,  Po = 208.9824e-3,  At = 209.9871e-3,  Rn = 222.0176e-3,
                Fr = 223.0197e-3, Ra = 226.0254e-3, Ac = 227.0279e-3, Th = 232.0381e-3,
                Pa = 231.03588e-3, U = 238.0508e-3, Np = 237.0482e-3, Pu = 244.0482e-3
            )


"""
get_element_weight(el)

function to return the element weights
#   Usage:
    get_element_weight(el)
-   'el::String' : Element name
"""
get_element_weight(el::String) = elementWeight[Symbol(titlecase(strip(el)))]            


"""
Function to create thermo object. The function reads the therm.dat file and 
parses the content based on the ideal gas object to create the thermo data object
The function returns SpeciesThermoObj
#   Usage:    
    create_thermo(species::Array{T,1}, thermo_file::AbstractString )
-   species  : Array of species names 
-   thermo_file : name of the thermo file including the path 
"""
function create_thermo(species::Array{T,1}, thermo_file::AbstractString ) where T <: AbstractString
    #Create an array of thermo objects which stores the thermo data of all species
    species_thermo = Dict{String,NASAThermo}()
    molecular_weights_vector = Array{Float64,1}()    
        
    open(thermo_file) do io
        lno = 0
        local  species_name, phase, mol_wt,LT,HT,CT
        comp = Dict{String,Int64}()
        local poly_coeff = Array{Float64,1}()
        while !eof(io)
            data_string = readline(io)     
            if length(data_string) != 0
                if SubString(data_string,1,3) == "END"
                    break
                else                                        
                    sp = uppercase(strip(data_string[1:min(length(data_string),18)]))
                    if sp in species
                        lno = 1
                        species_name, phase,comp, mol_wt,LT,HT,CT = parase_thermo_species_data(data_string)                        
                        lno += 1
                        empty!(poly_coeff)
                    elseif lno > 1 && lno < 5
                        append!(poly_coeff,parse_thermo_polynomials(lno,data_string))
                        lno += 1
                        if lno > 4
                            nasa_thermo = NASAThermo(species_name,phase,comp ,mol_wt,LT,HT,CT,poly_coeff[8:14],poly_coeff[1:7])                            
                            species_thermo[nasa_thermo.name] = nasa_thermo                            
                        end
                    end
                end
            end
        end
    end 
    species_thermo_array = Array{NASAThermo,1}()
    for sp in species
        push!(species_thermo_array,species_thermo[sp])
    end
    for td in species_thermo_array
        push!(molecular_weights_vector,td.molWt)
    end
    
    return SpeciesThermoObj(molecular_weights_vector, species_thermo_array)    
end

"""
Function for parsing the first line of thermo data.
    It returns the species name, phase, molecular weight, high 
    temperature, low temperature and Common temperature limits
    Not for external calls
#   Usage:
    parase_thermo_species_data(data_string) 
"""
function parase_thermo_species_data(data_string::AbstractString)
    species_name = uppercase(strip(SubString(data_string,1,18)))
    a1 = split(strip(SubString(data_string,25,29)))
    a2 = split(strip(SubString(data_string,30,34)))
    a3 = split(strip(SubString(data_string,35,39)))
    a4 = split(strip(SubString(data_string,40,44)))
    elements_data = (a1,a2,a3,a4)
    mol_wt::Float64 = 0
    composition = Dict{String,Int64}()
    for el in elements_data
        if !isempty(el)
            mol_wt += parse(Float64,el[2])*get_element_weight(String(el[1]))            
            composition[String(el[1])] = parse(Int64,el[2])
        end
    end
    phase = String(SubString(data_string,45,45))
    if length(phase)==0
        throw(error("Phase not specified for species $spname"))
    end
    LTs = String(SubString(data_string,46,57))
    if length(LTs) == 0
        throw(error("Low temperature limit not specified for species $spname"))
    end
    HTs = String(SubString(data_string,58,67))
    if length(HTs) == 0
        throw(error("High temperature limit not specified for species $spname"))
    end
    CTs = String(SubString(data_string,68,75))
    if length(CTs)==0
        throw(error("Common temperature limit not specified for species $spname"))
    end
    LT = parse(Float64, LTs)
    HT = parse(Float64, HTs)
    CT = parse(Float64, CTs)    
    return (species_name,phase, composition ,mol_wt,LT,HT,CT)
end


"""
This function will parse the polynomials
    Not for external calls
"""
function parse_thermo_polynomials(lno::Int64, thdata::AbstractString)
    poly_data = Array{Float64,1}()
    if lno == 2
        push!(poly_data,parse(Float64,String(SubString(thdata,1,15))))
        push!(poly_data,parse(Float64,String(SubString(thdata,16,30))))
        push!(poly_data,parse(Float64,String(SubString(thdata,31,45))))
        push!(poly_data,parse(Float64,String(SubString(thdata,46,60))))
        push!(poly_data,parse(Float64,String(SubString(thdata,61,75))))
    elseif lno == 3
        push!(poly_data,parse(Float64,String(SubString(thdata,1,15))))
        push!(poly_data,parse(Float64,String(SubString(thdata,16,30))))
        push!(poly_data,parse(Float64,String(SubString(thdata,31,45))))
        push!(poly_data,parse(Float64,String(SubString(thdata,46,60))))
        push!(poly_data,parse(Float64,String(SubString(thdata,61,75))))
    elseif lno == 4
        push!(poly_data,parse(Float64,String(SubString(thdata,1,15))))
        push!(poly_data,parse(Float64,String(SubString(thdata,16,30))))
        push!(poly_data,parse(Float64,String(SubString(thdata,31,45))))
        push!(poly_data,parse(Float64,String(SubString(thdata,46,60))))
    end
    return poly_data
end


"""
Calculate the specific heat of pure species J/mol-K
#   Usage-1:
    cp(thermo,T) 
-   thermo::NASAThermo: NASAThermo of the species
-   T::Float64: Temperature in K at which the property is required    
# Usage-2:
    cp(sp,T,thermo,ig)
-   sp::String : species name
-   T::Float64 : Temperature K
-   thermoObj::SpeciesThermoObj : Structure of SpeciesThermoObj
-   species_list::Array{String,1}  : List of species 
"""
function cp(thermo::NASAThermo, T::Float64)
    TVec = [T^i for i in 0:4]    
    T < thermo.cmt ? sum(thermo.ltp[1:5] .* TVec)*R : sum(thermo.htp[1:5] .* TVec)*R
end
cp(sp::String,T::Float64,thermoObj::SpeciesThermoObj,species_list::Array{String,1}) = cp(thermoObj.thermo_all[get_index(sp,species_list)],T)



"""
H(thermo::NASAThermo, T::Float64)
Calculates the enthalpy of pure species J/mol
#   Usage-1:
    H(thermo,T)
-   'thermo::NASAThermo': NASAThermo of the species
-   'T::Float64': Temperature in K at which the property is required   
# Usage-2:
    H(sp,T,thermo,ig)
-   sp::String : species name
-   T::Float64 : Temperature K
-   thermoObj::SpeciesThermoObj : Structure of SpeciesThermoObj
-   species_list::Array{String,1}  : List of species 
"""
function H(thermo::NASAThermo, T::Float64)
    TVec = [1, T/2.0, T^2/3.0, T^3/4.0, T^4/5.0, 1.0/T]    
    T < thermo.cmt ? sum(thermo.ltp[1:6] .* TVec)*R*T : sum(thermo.htp[1:6] .* TVec)*R*T
end
H(sp::String,T::Float64,thermoObj::SpeciesThermoObj,species_list::Array{String,1}) = H(thermoObj.thermo_all[get_index(sp,species_list)],T)


"""
S(thermo::NASAThermo, T::Float64)
Calculates the entropy of pure species J/mol-K
#   Usage-1:
    S(thermo,T)
-   thermo::NASAThermo: NASAThermo of the species
-   T::Float64: Temperature in K at which the property is required    
# Usage-2:
    S(sp,T,thermo,ig)
-   sp::String : species name
-   T::Float64 : Temperature K
-   thermo::SpeciesThermoObj : Structure of SpeciesThermoObj
-   species_list::Array{String,1}  : List of species 
"""
function S(thermo::NASAThermo, T::Float64)
    TVec = [log(T), T, T^2/2.0, T^3/3.0, T^4/4.0]    
    if T < thermo.cmt
        return (sum(thermo.ltp[1:5] .* TVec) + thermo.ltp[7])*R
    else
        return (sum(thermo.htp[1:5] .* TVec) + thermo.htp[7])*R
    end
end
S(sp::String,T::Float64,thermoObj::SpeciesThermoObj,species_list::Array{String,1}) = S(thermoObj.thermo_all[get_index(sp,species_list)],T)


"""
Calculates the specific heat of all species in J/mol-K
#   Usage
    cp_all(td,T)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K at which the property is rquired
"""
function cp_all(thermoObj::SpeciesThermoObj,T::Float64)
    return map(x->cp(x,T),thermoObj.thermo_all)
end


"""
Calculates the enthalpy of all species in J/mol
#   Usage
    H_all(td,T)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K at which the property is rquired
"""
function H_all(thermoObj::SpeciesThermoObj,T::Float64)
    return map(x->H(x,T),thermoObj.thermo_all)
end


"""
Calculates the entropy of all species in J/mol-K
#   Usage
    S_all(td,T)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K at which the property is rquired
"""
function S_all(thermoObj::SpeciesThermoObj,T::Float64)
    return map(x->S(x,T),thermoObj.thermo_all)
end

"""
Hmix(thermoObj::SpeciesThermoObj,T::Float64,mlf::Array{Float64,1})    
Calculates the enthalpy of a mixture J/mol
#   Usage
    Hmix(td,T,mlf)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K
-   'mlf::Array{Float64,1}' : species mole fractions    
"""
function Hmix(thermoObj::SpeciesThermoObj,T::Float64,mlf::Array{Float64,1})    
    all_species_h = H_all(thermoObj,T)        
    return sum(all_species_h .* mlf)    
end

"""
Calculates the specific heat of a mixture in J/mol-K
#   Usage
    cpmix(td,T,mlf)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K
-   'mlf::Array{Float64,1}' : species mole fractions    
"""
function cpmix(thermoObj::SpeciesThermoObj,T::Float64,mlf::Array{Float64,1})
    all_species_cp = S_all(thermoObj,T)
    return sum(all_species_cp .* mlf)
end

"""
Smix(thermoObj::SpeciesThermoObj,T::Float64,p::Float64,mlf::Array{Float64,1})
Calculates the entropy of a muxture in J/mol-K
#   Usage
    Smix(thermoObj,T,p,mlf)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K 
-   'p::Float64' : total pressure Pa
-   'mlf::Array{Float64,1} ': mole fractions
"""
function Smix(thermoObj::SpeciesThermoObj,T::Float64,p::Float64,mlf::Array{Float64,1})
    all_species_s = S_all(thermoObj,T)
    xp = mlf .* (p/p_std)
    #logterm = [x <= 0 ? 0 : log(x) for x in xp]
    logterm = collect(map(x->x <= 0 ? 0 : log(x) ,xp))
    return sum((all_species_s - (logterm .* R)) .* mlf)
end

"""
Calculates the Gibbs free energy of a muxture in J/mol
#   Usage
    Gmix(thermoObj,T,p,mlf)
-   'thermoObj::SpeciesThermoObj' : Structure of SpeciesThermoObj
-   'T::Float64' : Temperature in K 
-   'p::Float64' : total pressure Pa
-   'mlf::Array{Float64,1} ': mole fractions
"""
function  Gmix(thermoObj::SpeciesThermoObj, T::Float64,p::Float64,mlf::Array{Float64,1})
    hmix = Hmix(thermoObj,T,mlf)
    smix = Smix(thermoObj,T,p,mlf)
    return hmix - T*smix
end


#end of module IdealGas
end