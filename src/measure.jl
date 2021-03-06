typealias ID Integer

# check if the bit at pos is 1 or not
# return true if it is 1
function isone(bit::Integer,pos::Integer)
    return Bool((bit>>(pos-1))&1)
end

# return the value of given bits with ids in the given bit squence
# eg. |10101> id:(1,3) returns 0b11 = 3
function bitvalue(bit::Integer,id::ID...)
    res = 0
    n = 0
    for i in id
        if isone(bit,i)
            res+=1<<n
        end
        n+=1
    end
    return res
end

function flip(bit::Integer,pos::Integer)
    return bit$(1<<(pos-1))
end

function lowerbit(bit::Integer,id::ID...)
    for j in id
        isone(bit,j)?bit=flip(bit,j):nothing
    end
    return bit
end

function expandbit(bit::Integer,id::ID...)
    res = 0
    len = length(id)
    for i=1:len
        res += ((bit>>(i-1))&1)<<(id[i]-1)
    end

    return res
end

function measure(state::AbstractVector,gate::Gate,id::ID...)
    #bounds check
    @assert length(state) ≥ 2^length(id) "Out of Bounds"
    @assert length(id) == gate.bitnum "number of operational bits does not match gate"
    @assert issorted([id...]) "id should be input from lower index to higher index"
    #return immediatly if it's a single bit operation
    if length(state)==2
        return state|>gate
    end

    state_len = length(state)
    res       = zeros(Complex,state_len)
    
    for i=0:state_len-1# loop from |00..0> to |11..1>
        ground = lowerbit(i,id...)
        local_state_catch = state[i+1]*([x==bitvalue(i,id...)?1:0 for x=0:2^gate.bitnum-1]|>gate)

        state_catch = zeros(Complex,state_len)
        for j = 1:length(local_state_catch)
            state_catch[ground+expandbit(j-1,id...)+1] = local_state_catch[j]
        end
        res += state_catch
    end

    return res
end

function measure(
    state::AbstractVector,
    control::AbstractVector,
    gate::Gate,
    id::ID...)
    # consts
    state_len = length(state)
    numbits = Int(log2(state_len))

    res = zeros(Complex,state_len)
    for i = 0:state_len-1
        #AND
        flag = true
        for j in control
            flag = flag&&isone(i,j)
        end

        if flag
        # control bits are all true
            ground = lowerbit(i,id...)
            local_state_catch = state[i+1]*([x==bitvalue(i,id...)?1:0 for x=0:2^gate.bitnum-1]|>gate)

            state_catch = zeros(Complex,state_len)
            for j = 1:length(local_state_catch)
                state_catch[ground+expandbit(j-1,id...)+1] = local_state_catch[j]
            end

            res += state_catch
        else
        # at least one of the control bits is false
            res += [x==(i+1)?Complex(state[i+1]):0+0im for x=1:state_len]
        end
    end
    return res
end

function measure(cc::Circuit,state::AbstractVector)
    for g in cc.gates
        state = measure(state,g...)
    end
    return normalize(state)
end

export measure

# state = measure(state,Hadamard,1)
# @show state
# @show normalize(measure(state,Hadamard,2))
# @show InitState(2)
# @show norm(measure([0,1/sqrt(2),1/sqrt(2),0],[2],Hadamard,1))
