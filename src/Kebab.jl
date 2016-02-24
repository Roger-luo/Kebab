module Kebab

import Base:bin

include("consts.jl")
include("truthtable.jl")
include("mathfunc.jl")

include("Adia/Adiabatic.jl")
include("gates.jl")


export CoolingAssit,AdiaSystem,TruthTable

end

using Kebab

H = AdiaSystem([TruthTable(0b1001,[1,2])],2,1e3)
bitnum = 2
CoolingAssit(H,bitnum)