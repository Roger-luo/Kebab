function check(gate::Gate,id...)
    if typeof(gate.op) <: AbstractMatrix
        @assert size(gate.op)==(2^length(id),2^length(id)) "Error: Operational Matrix is not in Correct Size"
    elseif typeof(gate.op) <: Function
        try
            gate.op(InitState(length(id)))
        catch
            error("Operator Is Not For $(length(id)) Qubits")
        end
    end
end

function addgate!(
    circuit::Circuit,
    gate::Gate,
    id::Integer...
    )
    check(gate,id...)
    push!(circuit.gates,(gate,id...))
    circuit.bit_num = maximum(id)
end

function addgate!(
    circuit::Circuit,
    control::AbstractVector,
    gate::Gate,
    id::Integer...
    )
    check(gate,id...)
    push!(circuit.gates,(control,gate,id...))
    circuit.bit_num = maximum(id)
end