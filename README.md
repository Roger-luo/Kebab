# Kebab

Kebab is a quantum computation framework for Julia.The structure is to realize quantum algorithm simulation under quantum circuit model. And add different gates to conbine different computation method.

# current functions
- **Quantum Circuit** quantum circuit is basic structure to use this framework.[have not been finished]

- **Adiabatic Computation** Adiabatic Computation Module is realized by add adiabatic evolution gates

- **Cooling** the cooling method comes from the [doi:10.103](http://www.nature.com/nphoton/journal/v8/n2/full/nphoton.2013.354.html)

# User Guide
To use this package, use `Pkg.clone(https://github.com/Roger-luo/AdiaRoll.jl.git)`

or download the source codes from github.

# One example for cooling-assist computation

```julia
using Kebab
using PyPlot

function CoolingAssit(
    Hs::AdiaSystem,
    bitnum::Real;
    dt = 1e-2
    )
    #init variables
    const evotime = Hs.maxtime
    state = convert(Array{Complex,1},[1/sqrt(2^bitnum) for i=1:2^bitnum])

    eigens = TimeEvoModule!(state,Hs,1/3;dt=dt)

    #cooling module
    gamma,t = CoolingPara(Hs)
    print("enter cooling module 1\n")
    CoolingModule!(state,Hs,gamma,t;n=5)
    print("finish cooling\n")

    eigens = TimeEvoModule!(state,Hs,1/3;dt=dt)

    #cooling module
    gamma,t = CoolingPara(Hs)
    print("enter cooling module 2\n")
    CoolingModule!(state,Hs,gamma,t;n=5)
    print("finish cooling\n")

    eigens = [eigens TimeEvoModule!(state,Hs,1/3;dt=dt)]
    
    success_prob = abs2([x==findmin(diag(Hs.HP))[2]?1:0 for x=1:2^bitnum])
    return eigens,state,success_prob[1]
end

H = AdiaSystem([TruthTable(0b1001,[1,2])],2,1e3)
eigen,state,p = CoolingAssit(H,2)

print("pass\n")

figure(1)
for i=1:size(eigen)[1]
    plot(real(eigen[i,:]))
end

show()
```