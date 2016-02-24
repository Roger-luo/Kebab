include("consts.jl")
include("Adiabatic.jl")

function Hadamard(qubit::AbstractVector)
    return hadamard*qubit
end

function gamma(qubit::AbstractVector,gamma::Real)
    return [1 0;0 -im*exp(im*gamma)]*qubit
end
