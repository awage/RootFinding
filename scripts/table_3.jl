using DrWatson
@quickactivate
using CodecZlib
# using CairoMakie
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))



function print_table_all()
    β_range = range(-1,1, step  = 0.5)
    res = 100; ε = 1.e-14;  max_it = 50; 
    force = true

    open("table1_dat.txt","w") do io
    for i in  1:14, k in 1:length(fam_list)
        prefix = string("steph_tan_f",i, "_g",k)
        println(string_list[i])
        print(io,"{\\footnotesize ", string_list[i], "}" )

        # g(z) = tanh(abs(z))*exp(im*angle(z))
        N = stephenson_map(func_list[i], fam_list[k])

        # Mean iterations
        m_it = _get_mean_it(N, res, ε, max_it; prefix, force) 
        print(io," & ",  round(m_it, digits =1))

        # Non converging points 
        nc = _get_mean_nc(N, res, ε, max_it; prefix, force = false) 
        print(io, " & ",  round(Int,100-100*nc)) # print convergence percentage

        # Computational time 
        t0_ref = _get_mean_t0(N, res, ε, max_it; prefix, force = false)
        @show t0_ref
        print(io, " & ",  round(t0_ref*1e6, digits =2))

        @show q = _get_q(N, res, ε, max_it) 
        print(io, " & ", round(q, digits = 2))

        println(io," \\\\")
    end
    end
end

print_table_all()



