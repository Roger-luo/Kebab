function realtimeop(
    state::AbstractVector,
    Hs::AdiaSystem,
    dt::Real
    )
    res = trotter(-im*dt*(1-Hs.location)*Hs.HB,-im*dt*Hs.location*Hs.HP,3)*state
    Hs.location += dt/Hs.maxtime
    return res
end

function realtimeop!(
    state::AbstractVector,
    Hs::AdiaSystem,
    dt::Real
    )
    copy!(state,trotter(-im*(1-Hs.location)*Hs.HB,-im*Hs.location*Hs.HP,3)*state)
    Hs.location += dt/Hs.maxtime
end
