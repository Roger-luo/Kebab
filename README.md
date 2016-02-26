# Kebab

Kebab is a quantum computation framework for Julia.Kebab aims to realize quantum algorithm simulation in different level(eg.hardware level/software level/...) under different computation model and provide useful support for experiments. 

# current functions
- **Quantum Circuit** quantum circuit is basic structure to use this framework.[have not been finished]
- **Quantum Fourier Transformation** QFT is based on quantum circuit

- **Adiabatic Computation** Adiabatic Computation Module is realized by add adiabatic evolution gates

- **Cooling** the cooling method comes from the [doi:10.103](http://www.nature.com/nphoton/journal/v8/n2/full/nphoton.2013.354.html)

# Simple User Guide
To use this package, use `Pkg.clone(https://github.com/Roger-luo/Kebab.git)`

or download the source codes from github.

```
TruthTable(truthvalue,bitID)
```

use `TruthTable` to construct the cost function,for example a cost function define below

<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Baligned%7D%0A%26h_c%28z_1%2Cz_2%2Cz_3%29%20%3D%201%5Cquad%20%5Ctext%7Bif%20there%20is%20a%201%20among%20%7D%20%28z_1%2Cz_2%2Cz_3%29%5C%5C%0A%26h_c%28z_1%2Cz_2%2Cz_3%29%20%3D%200%5Cquad%20%5Ctext%7Bif%20not%7D%0A%5Cend%7Baligned%7D&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0" align="center" border="0" alt="\begin{aligned}&h_c(z_1,z_2,z_3) = 1\quad \text{if there is a 1 among } (z_1,z_2,z_3)\\&h_c(z_1,z_2,z_3) = 0\quad \text{if not}\end{aligned}" width="417" height="43" />

then the truth table of this cost function is

| bits value | truth value |
|------------|-------------|
| 0   0   0  |           0 |
| 0   0   1  |           1 |
| 0   1   0  |           1 |
| 0   1   1  |           0 |
| 1   0   0  |           1 |
| 1   0   1  |           0 |
| 1   1   0  |           0 |
| 1   1   1  |           0 |

so the truth table construct as following,the `truthvalue` contains the truth value as a number
and the `bitID` contains the ID of bits in this truth table

```julia
TruthTable(0b01101000,[1,2,3])
```

for multi-thread calculation,use the `blas_set_num_threads` method

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