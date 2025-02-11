using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))

function compute_figure(N, x, ε, max_it)
    ds = FunIterator(N, x)
    n = @timed _get_iterations!(ds, ε, max_it)
    xf, _ = get_state(ds) 
    return n, xf
end

function print_table_all()
    ε = 1.e-14;  max_it = 60; 
    setprecision(BigFloat, 50; base = 10)

    open("table4_dat.txt","w") do io
    for i in 1:28
        print(io,L"{\footnotesize $f_{", i, L"}$}" )

        # Iterations
        xf_v = []
        for k in 1:length(g_list)
            N = stephenson_map(F_list[i], g_list[k], length(X0[i]))
            n, xf = compute_figure(N, X0[i], ε, max_it)
            @show xf
            push!(xf_v, xf)
            it = n.value
            print(io," & ", it)
        end
        
        #  Final point 
         for k in 1:length(g_list)
             # N = stephenson_map(F_list[i], g_list[k])
             # n, xf = compute_figure(N, X0[i], ε, max_it)
             print(io," & ")
             for x in xf_v[k]; print(io, round(Float64(x), digits =1), " "); end
         end

        println(io," \\\\")
    end

    end
end

print_table_all()



