function Hadamard(qubit::AbstractVector)
    return hadamard*qubit
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

# function Deutsch(qubit::AbstractVector,theta::Real)