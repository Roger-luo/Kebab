module Kebab

import Base:bin

include("consts.jl")
include("truthtable.jl")
include("mathfunc.jl")

include("Adia/Adiabatic.jl")
include("gates.jl")


export 
    AdiaSystem,#Adiabatic System Object
    TruthTable,#TruthTable Object
    TimeEvoModule!,#Time evolution gate
    Hadamard,
    gamma,
    #operators
    realtimeop,#real time operator
    realtimeop!,
    #cooling
    cooling!,#single cooling module
    CoolingModule!,#multi-times cooling module
    CoolingPara,#choose cooling parameter
    #math functions
    trotter,
    normalize!,
    ⊗,⊕
end