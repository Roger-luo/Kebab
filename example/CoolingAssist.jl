using Kebab
using PyPlot

blas_set_num_threads(4)

function CoolingAssit(
    Hs::AdiaSystem,
    bitnum::Real;
    dt = 1e-2
    )
    #init variables
    const evotime = Hs.maxtime
    state = convert(Array{Complex,1},[1/sqrt(2^bitnum) for i=1:2^bitnum])
    prob = 1

    eigens = TimeEvoModule!(state,Hs,1/3;dt=dt)

    #cooling module
    gamma,t = CoolingPara(Hs)
    print("enter cooling module 1\n")
    prob = CoolingModule!(state,Hs,gamma,t;n=5)
    print("finish cooling\n")

    eigens = TimeEvoModule!(state,Hs,1/3;dt=dt)

    #cooling module
    gamma,t = CoolingPara(Hs)
    print("enter cooling module 2\n")
    prob *= CoolingModule!(state,Hs,gamma,t;n=5)
    print("finish cooling\n")

    eigens = [eigens TimeEvoModule!(state,Hs,1/3;dt=dt)]
    
    success_prob = prob*abs2([x==findmin(diag(Hs.HP))[2]?1:0 for x=1:2^bitnum])
    return eigens,state,success_prob[1]
end

# H = AdiaSystem([TruthTable(0b10010111,[1,2,3]),
#     TruthTable(0b10010111,[2,3,4]),
#     TruthTable(0b10010111,[3,4,5]),
#     TruthTable(0b10010111,[5,6,7]),
#     TruthTable(0b10010111,[1,5,6]),
#     TruthTable(0b10010111,[2,3,6]), 
#     ],7,1e3)

H = AdiaSystem([TruthTable(0b10010111,[1,2,3])],3,1e3)

eigen,state,p = CoolingAssit(H,3)

print("pass\n")

figure(1)
for i=1:size(eigen)[1]
    plot(real(eigen[i,:]))
end

show()