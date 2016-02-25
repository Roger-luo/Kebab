function Hadamard(state::AbstractVector)
    @assert length(state)==2 "Hadamard gate is a single-qubit gate"
    return hadamard*state
end

function Pauli_X(state::AbstractVector)
    @assert length(state)==2 "Pauli_X gate is a single-qubit gate"
    return sigmax*state
end

function Pauli_Y(state::AbstractVector)
    @assert length(state)==2 "Pauli_Y gate is a single-qubit gate"
    return sigmay*state
end

function Pauli_Z(state::AbstractVector)
    @assert length(state)==2 "Pauli_Z gate is a single-qubit gate"
    return sigmaz*state
end

function R_k(state::AbstractVector,k::Int64)
    @assert length(state)==2 "R-k gate is a single-qubit gate"
    return [1 0;0 e^(2*π*im/(2^k))]*state
end
# function gamma(qubit::AbstractVector,gamma::Real)
#     return [1 0;0 -im*exp(im*gamma)]*qubit
# end

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
    res = zeros(bitnum2)
    for i = 1:length(state)
        if isone(i-1,cid)&&(!isone(i-1,ucid))
            temp = zeros(bitnum2)
            single = gate*[1,0]
            temp[i] = state[i]*single[1]
            temp[flip(i-1,ucid)+1] = state[i]*single[2]
        elseif isone(i-1,cid)&&isone(i-1,ucid)
            temp = zeros(bitnum2)
            single = gate*[0,1]
            temp[i] = state[i]*single[2]
            temp[flip(i-1,ucid)+1] = state[i]*single[1]
        else
            temp = zeros(bitnum2)
            temp[i] = state[i]
        end
        res+=temp
    end
    return res
end


function ControlGate(
    gate::Function,
    state::AbstractVector,
    CID::Int64,
    UCID::Int64
    )
    bitnum2 = length(state)
    @assert (bitnum2 & (bitnum2-1))==0 "Wrong qubits state"

    bitnum = Int(log2(bitnum2))
    cid = bitnum - CID +1
    ucid = bitnum - UCID +1
    res = zeros(bitnum2)
    for i = 1:length(state)
        if isone(i-1,cid)&&(!isone(i-1,ucid))
            temp = zeros(bitnum2)
            single = gate([1,0])
            temp[i] = state[i]*single[1]
            temp[flip(i-1,ucid)+1] = state[i]*single[2]
        elseif isone(i-1,cid)&&isone(i-1,ucid)
            temp = zeros(bitnum2)
            single = gate([0,1])
            temp[i] = state[i]*single[2]
            temp[flip(i-1,ucid)+1] = state[i]*single[1]
        else
            temp = zeros(bitnum2)
            temp[i] = state[i]
        end
        res+=temp
    end
    return res
end

# function Deutsch(qubit::AbstractVector,theta::Real)