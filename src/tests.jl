using Kebab
state = [0,0,0,1]
@show ControlGate(hadamard,state,1,2)
@show blkdiag(sparse([1 0;0 1]),sparse(hadamard))*state