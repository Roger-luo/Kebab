function boundscheck(Hs::AdiaSystem,gamma::Real,t::Real)
    energy = eig(full(Hamiltonian(Hs)))[1]

    Min = minimum(energy)
    Max = maximum(energy)

    print("minimum energy : $(Min)\n")
    print("maximum energy : $(Max)\n")
    @show Min*t-gamma
    @show Max*t-gamma
    @show (gamma-pi/2)/Min
    @show (gamma+pi/2)/Min
    @show (gamma-pi/2)/Max
    @show (gamma+pi/2)/Max
    @show t
end