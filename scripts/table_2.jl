using DrWatson
@quickactivate
using LaTeXStrings
using Statistics
using JLD2
using Roots
using LinearAlgebra:norm
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))


function beta(j,β; res = 500, ε = 1e-4, max_it = 50)
    N_β = beta_map_anneal(func_list[j]); 
    d = @dict(N_β, β, res, ε, max_it) # parametros
    d1 = compute_basins_prox(d)
    @unpack grid, basins, iterations, attractors = d1
    ind = findall(basins .>= 1)
    # Only take into account ICs that have converged
    m1 = mean(iterations[ind])
    f_div = length(ind)/res^2
    return m1, f_div
end

function compute_annealing()
    ε = 1e-15
    beta_it_0 = zeros(15)
    beta_it_1 = zeros(15)
    beta_it_an = zeros(15)
    beta_it_an2 = zeros(15)
    for n in 1:15
        # print(string_list[n], " & ")
        # for (k,β) in enumerate(beta_range)
            β0(x,y) = 0. 
            m2, f2 = beta(n,β0;ε); 
            beta_it_0[n] = m2
            # print(round(m2,digits =1) , " & ")
            β1(x,y) = 1.
            m2, f2 = beta(n,β1;ε); 
            beta_it_1[n] = m2
            # print(round(m2,digits =1) , " & ")
            β3(x,y) =  (1 - 1/(1+(0.1*abs(x))^6))
            m3, f3 = beta(n,β3;ε); 
            beta_it_an[n] = m3
            # print(round(m3,digits =1), " & ")
            β4(x,y) =  2*x/(x+y)
            m4, f4 = beta(n,β4;ε); 
            beta_it_an2[n] = m4
            # print(round(m4,digits =1))
        # end
        # println(" \\\\")
        # println("\\hline")
    end
    @save "newton_annealing.jld2" beta_it_0 beta_it_1 beta_it_an beta_it_an2
end


function print_table()
    @load  "newton_annealing.jld2" beta_it_0 beta_it_1 beta_it_an beta_it_an2
    for n in 1:15
        print(string_list[n], " & ")
        print(round(beta_it_0[n],digits =1) , " & ")
        print(round(beta_it_1[n],digits =1) , " & ")
        # print(round(beta_it_an[n],digits =1), " & ")
        print(round(beta_it_an2[n],digits =1))
        println(" \\\\")
        println("\\hline")
    end
end

# compute_annealing()
print_table()
