module Adia

import Base:bin

include("consts.jl")
include("truthtable.jl")
include("mathfunc.jl")

type AdiaSystem
    HB::AbstractMatrix
    HP::AbstractMatrix
    location::Real
    maxtime::Real
    bitnum::Int64

    function AdiaSystem(expr::LogicExpr,bitnum::Int,maxtime::Real)
        HB = bHamilton(bitnum)
        HP = pHamilton(expr,bitnum)
        location = 0
        
        new(HB,HP,location,maxtime,bitnum)
    end
end

include("AdiaHamilton.jl")

function realtimeop(
    state::AbstractVector,
    Hs::AdiaSystem;
    dt=1e-2
    )
    res = trotter(-im*dt*(1-Hs.location)*Hs.HB,-im*dt*Hs.location*Hs.HP,3)*state
    Hs.location += dt/Hs.maxtime
    return res
end

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

function cooler!(state::AbstractVector,Hs::AdiaSystem,gamma::Real,t::Real)
    # boundscheck(Hs,gamma,t)
    energy = eig(full(Hamiltonian(Hs)))[1]
    @assert -pi/2<minimum(energy)*t-gamma<pi/2 "bad cooling parameters"
    @assert -pi/2<maximum(energy)*t-gamma<pi/2 "bad cooling parameters"

    copy!(state,normalize!((state-im*exp(im*gamma)*
            trotter(-im*t*(1-Hs.location)*Hs.HB,
                -im*t*Hs.location*Hs.HP,3)*state)/2))
end

function heater!(state::AbstractVector,Hs::AdiaSystem,gamma::Real,t::Real)
    # boundscheck(Hs,gamma,t)
    energy = eig(full(Hamiltonian(Hs)))[1]
    @assert -pi/2<minimum(energy)*t-gamma<pi/2 "bad cooling parameters"
    @assert -pi/2<maximum(energy)*t-gamma<pi/2 "bad cooling parameters"

    copy!(state,normalize!((state+im*exp(im*gamma)*
                trotter(-im*t*(1-Hs.location)*Hs.HB,
                    -im*t*Hs.location*Hs.HP,3)*state)/2))
end


function realtimeop!(
    state::AbstractVector,
    Hs::AdiaSystem,
    dt
    )
    copy!(state,trotter(-im*(1-Hs.location)*Hs.HB,-im*Hs.location*Hs.HP,3)*state)
    Hs.location += dt/Hs.maxtime
end

function cooling!(
    state::AbstractVector,
    Hs::AdiaSystem,
    gamma::Real,
    t::Real
    )
    dice = rand()
    if dice <= 0.5*(1-sin(gamma))
        #get |1> (higher energy)
        heater!(state,Hs,gamma,t)
        return false
    else
        #get |0> (lower energy)
        cooler!(state,Hs,gamma,t)
        return true
    end
end

function CoolingModule!(
    state::AbstractVector,
    Hs::AdiaSystem,
    gamma::Real,
    t::Real;
    n=5)

    count = 0

    while count<n
        flag = cooling!(state,Hs,gamma,t)
        if flag
            count += 1
        end
    end
end

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

function CoolingPara(Hs::AdiaSystem)
    Eigen = eig(full(Hamiltonian(Hs)))[1]
    maxEigen = maximum(Eigen)
    minEigen = minimum(Eigen)

    gamma = (maxEigen+minEigen)/(maxEigen-minEigen) * pi/2 * 0.1

    if (gamma-pi/2)>0
        t = 0.5*((gamma-pi/2)/minEigen + (gamma+pi/2)/maxEigen)
    else
        t = 0.5*((gamma-pi/2)/maxEigen + (gamma+pi/2)/maxEigen)
    end
    return gamma,t
end

function evolution(
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

export AdiaSystem,TruthTable,realtimeop,evolution,CoolingModule!

end

using Adia

H = AdiaSystem([TruthTable(0b1001,[1,2])],2,1e3)

bitnum = 2

evolution(H,bitnum)

# state = [1/sqrt(2^bitnum) for i=1:2^bitnum]

# CoolingModule!(state,H,pi/2,1e-3;n=5)


# H = AdiaSystem([0.5 0;0 -0.5],[1/sqrt(2) 1/sqrt(2);1/sqrt(2) -1/sqrt(2)],0,1e3)

# s = [1/sqrt(2),1/sqrt(2)]

# for i=1:30
#     s = cooler(s,H,0,1.5)
# end
# @show s
# # @show realtimeop([1,0],H)
# # @show H