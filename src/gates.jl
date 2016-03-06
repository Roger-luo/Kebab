function R_k(k::Int64)
    return [1 0;0 e^(2*π*im/(2^k))]
end

# basic gates
#wrapper for some commonly used gates
Hadamard = Gate(hadamard,1)
PI_8    = Gate(R_k(4),1)

Pauli_X = Gate(sigmax,1)
Pauli_Y = Gate(sigmay,1)
Pauli_Z = Gate(sigmaz,1)