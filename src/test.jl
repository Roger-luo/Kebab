function normalize!(state::AbstractVector)
    copy!(state,state/norm(state))
end

H = [0.5 0;0 -0.5]

function cooler(state::AbstractVector,gamma::Real,t::Real)
    return normalize!(0.5*(I-im*exp(im*gamma)*expm(-im*t*H))*state)
end

function imagtime(state)
    return normalize!(expm(-H*1)*state)
end

s = [1/sqrt(2),1/sqrt(2)]

for i=1:1000
    s = cooler(s,0,pi/2)
    # s = imagtime(s)
end