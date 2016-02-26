function SingleQubitQFT!(state::AbstractVector,bitID::Integer)
    state_len = length(state)
    bit_num = Int(log2(state_len))

    @assert bitID<=bit_num "bit ID out of bounds"
    
    copy!(state,Hadamard(state,bitID))
    for i=1:(bitID-1)
        copy!(state,C_R_k(i))
    end
end

function QFT!(state::AbstractVector)
    state_len = length(state)
    bit_num = Int(log2(state_len))
    
    for i=bit_num:1
        SingleQubitQFT!(state,i)
    end
end