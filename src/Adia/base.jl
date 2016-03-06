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

export AdiaSystem