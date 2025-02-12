using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function compute_figure(N, x, ε, max_it)
    ds = FunIterator(N, x)
    n, yy = _get_iterations!(ds, ε, max_it)
    xf, _ = get_state(ds) 
    if 5 ≤ n < max_it
        q = estimate_ACOC!(ds, n, yy)
    else
        q = 0
    end
    return n, xf, q
end

function print_table_all()
    ε = 1e-25;  max_it = 100; 
    setprecision(BigFloat, 50; base = 10)

    open("table4_dat.txt","w") do io
    for i in 1:28
        print(io,L"{\footnotesize $f_{", i, L"}$}" )

        # Iterations
        xf_v = []
        q_v = []
        for k in 1:length(g_list)
            N = stephenson_map(F_list[i], g_list[k], length(X0[i]))
            n, xf, q = compute_figure(N, X0[i], ε, max_it)
            @show xf
            push!(xf_v, xf)
            push!(q_v, q)
            # it = n.value
            print(io," & ", n)
        end
        
        #  Final point 
         for k in 1:length(g_list)
             if length(xf_v[k]) > 1
                 print(io," & (")
             else
                 print(io," & ")
             end
             for x in xf_v[k]; print(io, round(Float64(x), digits =1), ", "); end
             if length(xf_v[k]) > 1
                 print(io,")")
             end
         end

         # convergence order. 
         for k in 1:length(g_list)
             print(io," & ", round(Float64(q_v[k]), digits =1))
         end

        println(io," \\\\")
    end

    end
end

print_table_all()



