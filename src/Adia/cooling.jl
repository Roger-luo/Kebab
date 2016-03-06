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

function cooling!(
    state::AbstractVector,
    Hs::AdiaSystem,
    gamma::Real,
    t::Real
    )
    dice = rand()
    if dice <= 0.5*(1+sin(gamma))
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
    prob = 1

    while count<n
        flag = cooling!(state,Hs,gamma,t)
        if flag
            prob *= 0.5*(1+sin(gamma))
        else
            prob *= 0.5*(1-sin(gamma))
        end
        count += 1
    end
    return prob
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

export CoolingPara,CoolingModule,CoolingModule!