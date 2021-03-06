type TruthTable
    value::Integer
    bitID::AbstractVector{Integer}
end

typealias LogicExpr AbstractVector{TruthTable}

# function ClauseGenerator()

#return a bin string of given length
function bin(int::Integer,n::Integer)
    @assert length(bin(int))<=n

    res = ""
    for i=n-length(bin(int))
        res="$(res)0"
    end
    return res
end

# check the value of an assignment of bits
function TruthTableValue(BitsValue::Integer,truth_table::Integer,tablelen::Int64)
    TruthValue = bin(truth_table,2^tablelen)
    # @show BitsValue
    # @show TruthValue
    return parse(Int,TruthValue[BitsValue+1])
end

#expand a truth table
function TruthTableExpand(truth_table::Integer,bitID::AbstractVector{Integer},bitnum::Integer)
    # @assert bitnum<length(bin(truth_table))
    @assert length(bitID)<=bitnum
    sort!(bitID)
    res = ""
    for i = 0:2^bitnum-1
        strbin = bin(i,bitnum)
        strtemp = "0b"
        for id in bitID
            strtemp = "$(strtemp)$(strbin[id])"
        end
        res = "$(res)$(TruthTableValue(parse(Int,strtemp),truth_table,length(bitID)))"
    end
    return res
end

"""
```
State2TruthTable(state::AbstractVector,bitnum::Integer)
```
transfer the qubit state vector to result

Parameters
---
- `state::AbstractVector` a state vector
- `bitnum::Integer` number of qubits
"""
function State2TruthTable(state::AbstractVector,bitnum::Integer)
    res = abs(real(state))
    res = convert(Array{Int64,1},round(res))
    count=0
    for i in res
        if i==1
            print("$(bin(count,bitnum))\n")
        end
        count+=1
    end
end

# """
# Exact Cover Truth Table Generator
# ---
# """
# function ECGenerator()

export TruthTable