using DrWatson
@quickactivate
using CodecZlib
using CairoMakie
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))


function print_table_all()
    β_range = range(-1,1, step  = 0.5)
    res = 1000; ε = 1.e-14;  max_it = 50

    open("table_convergence.txt","w") do io
    for i in  [1:13; 15; 17]
        println(string_list[i])
        print(io,"{\\footnotesize ", string_list[i], "}" )
        N_β = beta_map(func_list[i])

        q = [_get_q(N_β, β, i, res, ε, max_it) for β in β_range]
        for k in eachindex(β_range)
            print(io, " & ", round(q[k], digits = 2))
        end
        println(io," \\\\")
    end
    end
end

print_table_all()

