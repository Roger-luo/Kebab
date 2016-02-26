function SingleQubitGate(
    gate::AbstractMatrix,
    state::AbstractVector,
    ID::Integer
    )
    @assert size(gate)==(2,2) "Input gate is not a single-qubit gate"

    if length(state)==2
        return gate*state
    end

    state_len = length(state)
    res = zeros(Complex,state_len)
    id = Int(log2(state_len))-ID+1

    for i=1:state_len
        single_state_temp = zeros(Complex,state_len)
        if isone(i-1,id)
            # |1>
            single = gate*[0,1]
            single_state_temp[i] = state[i]*single[2]
            single_state_temp[flip(i-1,id)+1] = state[i]*single[1]
        else
            # |0>
            single = gate*[1,0]
            single_state_temp[i] = state[i]*single[1]
            single_state_temp[flip(i-1,id)+1] = state[i]*single[2]
        end
        res+=single_state_temp
    end
    return res
end

function ControlGate(
    gate::AbstractMatrix,
    state::AbstractVector,
    CID::Int64,
    UCID::Int64
    )
    @assert size(gate) == (2,2) "The applied gate is not a single-bit gate"
    bitnum2 = length(state)
    @assert (bitnum2 & (bitnum2-1))==0 "Wrong qubits state"

    bitnum = Int(log2(bitnum2))
    cid = bitnum - CID +1
    ucid = bitnum - UCID +1
    res = zeros(Complex,bitnum2)
    for i = 1:length(state)
        if isone(i-1,cid)&&(!isone(i-1,ucid))
            temp = zeros(Complex,bitnum2)
            single = gate*[1,0]
            temp[i] = state[i]*single[1]
            temp[flip(i-1,ucid)+1] = state[i]*single[2]
        elseif isone(i-1,cid)&&isone(i-1,ucid)
            temp = zeros(Complex,bitnum2)
            single = gate*[0,1]
            temp[i] = state[i]*single[2]
            temp[flip(i-1,ucid)+1] = state[i]*single[1]
        else
            temp = zeros(Complex,bitnum2)
            temp[i] = state[i]
        end
        res+=temp
    end
    return res
end

## May not need Function type
# function ControlGate(
#     gate::Function,
#     state::AbstractVector,
#     CID::Int64,
#     UCID::Int64
#     )
#     bitnum2 = length(state)
#     @assert (bitnum2 & (bitnum2-1))==0 "Wrong qubits state"

#     bitnum = Int(log2(bitnum2))
#     cid = bitnum - CID +1
#     ucid = bitnum - UCID +1
#     res = zeros(bitnum2)
#     for i = 1:length(state)
#         if isone(i-1,cid)&&(!isone(i-1,ucid))
#             temp = zeros(bitnum2)
#             single = gate([1,0])
#             temp[i] = state[i]*single[1]
#             temp[flip(i-1,ucid)+1] = state[i]*single[2]
#         elseif isone(i-1,cid)&&isone(i-1,ucid)
#             temp = zeros(bitnum2)
#             single = gate([0,1])
#             temp[i] = state[i]*single[2]
#             temp[flip(i-1,ucid)+1] = state[i]*single[1]
#         else
#             temp = zeros(bitnum2)
#             temp[i] = state[i]
#         end
#         res+=temp
#     end
#     return res
# end

# basic gates
function R_k(k::Int64)
    return [1 0;0 e^(2*π*im/(2^k))]
end

#wrapper for some commonly used gates
function Hadamard(state::AbstractVector,ID::Integer)
    return SingleQubitGate(hadamard,state,ID)
end

function Pauli_X(state::AbstractVector,ID::Integer)
    return SingleQubitGate(sigmax,state,ID)
end

function Pauli_Y(state::AbstractVector,ID::Integer)
    return SingleQubitGate(sigmay,state,ID)
end

function Pauli_Z(state::AbstractVector,ID::Integer)
    return SingleQubitGate(sigmaz,state,ID)
end

#Control Gates

function C_R_k(
    state::AbstractVector,
    CID::Integer,UCID::Integer,
    k::Integer
    )
    return ControlGate(R_k(k),state,CID,UCID)
end
# function Deutsch(qubit::AbstractVector,theta::Real)