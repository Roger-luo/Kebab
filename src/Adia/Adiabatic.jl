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
    @assert 0<=Hs.location+evopercentage <=1 "evolutoin percentage out of bounds(should be in [0,1])"

    const evotime = Hs.maxtime
    for i=Hs.location*evotime:dt:(Hs.location+evopercentage)*evotime
        realtimeop!(state,Hs,dt)
    end
end

function CoolingAssit(
    Hs::AdiaSystem,
    bitnum::Real;
    dt = 1e-2
    )
    #init variables
    const evotime = Hs.maxtime
    state = convert(Array{Complex,1},[1/sqrt(2^bitnum) for i=1:2^bitnum])

    #data buffer
    eigentemp = eig(full(Hamiltonian(0,Hs.HB,Hs.HP)))
    eigenvalue = eigentemp[1]
    eigennum = length(eigenvalue)

    #start evolution

    @show Hs.location

    for i=0:dt:evotime/3
        realtimeop!(state,Hs,dt)

        eigentemp = eig(full(Hamiltonian(Hs)))
        append!(eigenvalue,eigentemp[1])
    end

    #cooling module
    gamma,t = CoolingPara(Hs)
    print("enter cooling module 1\n")
    CoolingModule!(state,Hs,gamma,t;n=5)

    print("finish cooling\n")

    for i=evotime/3:dt:2*evotime/3
        realtimeop!(state,Hs,dt)

        eigentemp = eig(full(Hamiltonian(Hs)))
        append!(eigenvalue,eigentemp[1])
    end

    #cooling module
    gamma,t = CoolingPara(Hs)
    print("enter cooling module 2\n")
    CoolingModule!(state,Hs,gamma,t;n=5)

    print("finish cooling\n")

    for i=2*evotime/3:dt:evotime
        realtimeop!(state,Hs,dt)

        eigentemp = eig(full(Hamiltonian(Hs)))
        append!(eigenvalue,eigentemp[1])
    end

    success_prob = abs2([x==findmin(diag(Hs.HP))[2]?1:0 for x=1:2^bitnum])
    return reshape(eigenvalue,eigennum,Int(length(eigenvalue)/eigennum)),state,success_prob[1]
end