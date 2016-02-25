type AdiaSystem
    HB::AbstractMatrix
    HP::AbstractMatrix
    location::Real
    maxtime::Real
    bitnum::Int64

    function AdiaSystem(expr::LogicExpr,bitnum::Int,maxtime::Real)
        HB = bHamilton(bitnum)
        HP = pHamilton(expr,bitnum)
        location = 0
        
        new(HB,HP,location,maxtime,bitnum)
    end
end

type Qubit
    state::AbstractVector
    bitID::Int64

    function Qubit(state::AbstractVector,bitID::Int64)
        @assert length(state) == 2 "state of qubit should be a 2 variable array"
        new(state,bitID)
    end
end