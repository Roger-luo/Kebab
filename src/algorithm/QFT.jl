function SingleQubitQFT!(state::AbstractVector,bitID::Integer)
    state_len = length(state)
    bit_num = Int(log2(state_len))

    @assert bitID<=bit_num "bit ID out of bounds"
    
    copy!(state,Hadamard(state,bitID))
    for i=1:(bitID-1)
        copy!(state,C_R_k(state,bitID-i,bitID,i))
    end
    return state
end

function SingleQubitQFT(state::AbstractVector,bitID::Integer)
    state_len = length(state)
    bit_num = Int(log2(state_len))

    @assert bitID<=bit_num "bit ID out of bounds"

    res = Hadamard(state,bitID)
    for i=1:(bitID-1)
        res = C_R_k(res,bitID-i,bitID,i)
    end
    return res
end

function QFT!(state::AbstractVector)
    state_len = length(state)
    bit_num = Int(log2(state_len))

    i = bit_num
    while i>0
        SingleQubitQFT!(state,i)
        i-=1
    end
    return state
end

function QFT(state::AbstractVector)
    state_len = length(state)
    bit_num = Int(log2(state_len))

    res = state
    i = bit_num
    while i>0
        res = SingleQubitQFT(res,i)
        i-=1
    end
    return res
end