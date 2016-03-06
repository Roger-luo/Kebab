module Kebab

# __precompile__()

import Base:bin

include("consts.jl")
include("base.jl")
include("construct.jl")
include("measure.jl")
include("gates.jl")

# include("consts.jl")
# include("mathfunc.jl")

include("Adia/truthtable.jl")
include("Adia/Adiabatic.jl")
# include("gates.jl")
# include("algorithm/QFT.jl")


# export 
#     #consts
#     σ₀,σ₁,σ₂,σ₃,hadamard,
#     #AdiaComput
#     AdiaSystem,#Adiabatic System Object
#     TruthTable,#TruthTable Object
#     TimeEvoModule!,#Time evolution gate
#     #circuit
#     Hadamard,R_k,C_R_k,
#     Pauli_X,Pauli_Y,Pauli_Z,
#     ControlGate,
#     SingleQubitGate,
#     #algorithm
#     QFT!,QFT,
#     #operators
#     realtimeop,#real time operator
#     realtimeop!,
#     #cooling
#     cooling!,#single cooling module
#     CoolingModule!,#multi-times cooling module
#     CoolingPara,#choose cooling parameter
#     #math functions
#     trotter,
#     ⊗,⊕
end

# using Kebab

# test = Circuit()
# addgate!(test,Hadamard,1)
# addgate!(test,Hadamard,2)
# state = InitState(2)
# @show state
# @show measure(test,state)