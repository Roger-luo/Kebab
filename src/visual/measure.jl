# check if the bit at pos is 1 or not
# return true if it is 1
function isone(bit::Integer,pos::Integer)
    return Bool((bit>>(pos-1))&1)
end

function bitvalue(bit::Integer,id::Integer...)
    res = 0
    for i in id
        if isone(bit,i)
            res+=1<<(i-1)
        end
    end
    return res
end

function flip(bit::Integer,pos::Integer)
    return bit$(1<<(pos-1))
end

function lowerbit(bit::Integer,id::Integer...)
    for j in id
        isone(bit,j)?bit=flip(bit,j):nothing
    end
    return bit
end

function expandbit(bit::Integer,id...)
    res = 0
    len = length(id)
    for i=1:len
        res += ((bit>>(i-1))&1)<<(id[i]-1)
    end

    return res
end

function measure(state::AbstractVector,gate::Gate,id::Integer...)
    #bounds check
    @assert length(state) ≥ 2^length(id) "Out of Bounds"
    @assert length(id) == gate.bitnum "number of operational bits does not match gate"
    @assert issorted([id...]) "id should be input from lower index to higher index"
    #return immediatly if it's a single bit operation
    if length(state)==2
        return state|>gate
    end

    state_len = length(state)
    bits      = log2(length(state))
    res       = zeros(Complex,state_len)
    
    for i=0:state_len-1# loop from |00..0> to |11..1>
        ground = lowerbit(i,id...)
        local_state_catch = state[i]*[x==bitvalue(i,id...)?1:0 for x=0:2^gate.bitnum-1]|>gate
        state_catch = zeros(Complex,state_len)
        for j = 1:length(local_state_catch)
            state_catch[ground+expandbit(j-1)] = local_state_catch[j]
        end
        res += state_catch
    end

    return res
end


# function measure(cc::Circuit)
