function R_k(k::Int64)
    return [1 0;0 e^(2*π*im/(2^k))]
end

# basic gates
#wrapper for some commonly used gates
Hadamard = Gate("H",hadamard,1)
PI_8    = Gate("π/8",R_k(4),1)
R_1 = Gate("R_1",R_k(1),1)
R_2 = Gate("R_2",R_k(2),1)
R_3 = Gate("R_3",R_k(3),1)
R_4 = Gate("R_4",R_k(4),1)
R_5 = Gate("R_5",R_k(5),1)
R_6 = Gate("R_6",R_k(6),1)

Pauli_X = Gate("σx",sigmax,1)
Pauli_Y = Gate("σy",sigmay,1)
Pauli_Z = Gate("σz",sigmaz,1)

export R_k,Hadamard,Pauli_Z,Pauli_Y,Pauli_X