include("base.jl")
include("Hamilton.jl")
include("operator.jl")
include("cooling.jl")

function TimeEvoModule!(
    state::AbstractVector,
    Hs::AdiaSystem,
    evopercentage::Real;
    dt=1e-2
    )
    @assert (0<=Hs.location+evopercentage)&&( Hs.location+evopercentage-1<1e-4) "evolutoin percentage out of bounds(should be in [0,1])"

    const evotime = Hs.maxtime
    eigens = eig(full(Hamiltonian(Hs)))[1]
    for i=Hs.location*evotime:dt:(Hs.location+evopercentage)*evotime
        realtimeop!(state,Hs,dt)

        eigens = [eigens eig(full(Hamiltonian(Hs)))[1]]
    end
    return eigens
end

export TimeEvoModule!