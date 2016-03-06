# Basic Definiations for Quantum Circuits
# Abstract Types
abstract AbstractGate
abstract AbstractOperator
abstract AbstractCircuit

typealias QuBit Array{Int,1}


# math functions
function isunitary(umatrix::AbstractMatrix)
    if inv(umatrix)==umatrix.'
        return true
    else
        return false
    end
end

type UnitaryMatrix <: AbstractOperator
    op::AbstractMatrix

    function UnitaryMatrix(umatrix::AbstractMatrix)
        if isunitary(umatrix)
            new(umatrix)
        else
            print("It is not a Unitary Matrix\n")
        end
    end
end

#trotter expansion
function trotter(A::AbstractMatrix,B::AbstractMatrix,P::Int64)
    return (expm(full(A/(2*P)))*expm(full(B/P))*expm(full(A/(2*P))))^P
end

function (⊗)(A::AbstractMatrix,B::AbstractMatrix)
    @assert size(A)[1]==size(A)[2]
    @assert size(B)[1]==size(B)[2]

    return kron(A,B)
end

function (⊕)(A::AbstractMatrix,B::AbstractMatrix)
    return full(blkdiag(sparse(A),sparse(B)))
end

##################################

function InitState(n::Integer=2)
    return [1/sqrt(2^n) for i=1:2^n]
end


immutable Gate{T}<:AbstractGate
    op::T
    bitnum::Integer

    Gate(op::T,n::Integer)=new(op,n)
end

Gate(op::AbstractMatrix,n::Integer)=Gate{AbstractMatrix}(op,n)
Gate(op::UnitaryMatrix,n::Integer)=Gate{UnitaryMatrix}(op,n)
Gate(op::Function,n::Integer)=Gate{Function}(op,n)

import Base.|>

function (|>)(state::AbstractVector,gate::Gate)
    #check bounds
    @assert 2^gate.bitnum == length(state) "This Gate is a $(gate.bitnum)-bit gate"

    if typeof(gate.op)<:Function
        return gate.op(state)
    elseif typeof(gate.op)<:AbstractMatrix
        return gate.op*state
    else
        error("Wrong Type For Quantum Gate")
    end
end

function (→)(state::AbstractVector,gate::Gate)
    #check bounds
    @assert 2^gate.bitnum == length(state) "This Gate is a $(gate.bitnum)-bit gate"

    if typeof(gate.op)<:Function
        return gate.op(state)
    elseif typeof(gate.op)<:AbstractMatrix
        return gate.op*state
    else
        error("Wrong Type For Quantum Gate")
    end
end

typealias Gates Array{Tuple,1}

type Circuit<:AbstractCircuit
    gates::Gates
    bit_num::Integer

    Circuit(gates::Gates,bit_num::Integer)=new(gates,bit_num)
end

Circuit()=Circuit(Array(Tuple,0),0)
Circuit(bit_num::Integer)=Circuit(Tuple,bit_num)