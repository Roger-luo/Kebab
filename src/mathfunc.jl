#trotter expansion
function trotter(A::AbstractMatrix,B::AbstractMatrix,P::Int64)
    return (expm(full(A/(2*P)))*expm(full(B/P))*expm(full(A/(2*P))))^P
end

# function normalize!(state::AbstractVector)
#     copy!(state,state/norm(state))
# end

function (⊗)(A::AbstractMatrix,B::AbstractMatrix)
    @assert size(A)[1]==size(A)[2]
    @assert size(B)[1]==size(B)[2]

    return kron(A,B)
end

function (⊕)(A::AbstractMatrix,B::AbstractMatrix)
    return full(blkdiag(sparse(A),sparse(B)))
end

function flip(bit::Integer,pos::Integer)
    return bit$(1<<(pos-1))
end

function isone(bit::Integer,pos::Integer)
    return Bool((bit>>(pos-1))&1)
end