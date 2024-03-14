using DrWatson
@quickactivate
using LaTeXStrings
using Statistics
using JLD2
include("../src/basins_compute.jl")
include("../src/function_stuff.jl")

function beta(j,β; res = 300, ε = 1e-4)
    N_β = beta_map_df(func_list[j]); 
    d = @dict(N_β, β, res, ε) # parametros
    d1 = compute_basins_prox(d)
    @unpack grid, basins, iterations, attractors = d1
    ind = findall(basins .>= 1)
    # Only take into account ICs that have converged
    m1 = mean(iterations[ind])
    f_div = length(ind)/res^2
    return m1, f_div
end

function compute_annealing()
    ε = 1e-10
    beta_it_0 = zeros(15)
    beta_it_1 = zeros(15)
    beta_it_an = zeros(15)
    for n in 1:15
        print(string_list[n], " & ")
        # for (k,β) in enumerate(beta_range)
            β0(x,y) = 0. 
            m2, f2 = beta(n,β0;ε); 
            beta_it_0[n] = m2
            print(round(m2,digits =1) , " & ")
            β1(x,y) = 1.
            m2, f2 = beta(n,β1;ε); 
            beta_it_1[n] = m2
            print(round(m2,digits =1) , " & ")
            β3(x,y) =  (1 - 1/(1+(0.1*abs(x))^6))
            # β3(x,y) =  abs(x) < 1e-2 ? 0 : 1
            m3, f3 = beta(n,β3;ε); 
            beta_it_an[n] = m3
            print(round(m3,digits =1))
        # end
        println(" \\\\")
        println("\\hline")
    end
    @save "newton_annealing.jld2" beta_it_0 beta_it_1 beta_it_an
end


function estimate_f14()
    T = 1500
    Nsim = 300^2
    β0(x) = 0.; 
    β1(x) = 1.; 
    βan(x) = (1 - 1/(1+(0.1*abs(x))^6))
    @show m0,v0 = estimate_Nit_real(14, β0, Nsim, T, 1e-10)
    @show m1,v1 = estimate_Nit_real(14, β1, Nsim, T, 1e-10)
    @show m1,v1 = estimate_Nit_real(14, βan, Nsim, T, 1e-10)
end

