using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function compute_figure(ds, ε, max_it)
    n, yy = _get_iterations!(ds, ε, max_it)
    xf, fx = get_state(ds) 
    @show xf, fx
    if 5 ≤ n < max_it
        q = estimate_ACOC!(n, yy)
    else
        q = 0
    end
    return n, xf, q
end

function print_table_all()
    ε = 1e-8;  max_it = 1000; 
    setprecision(BigFloat, 50; base = 10)

    open("table7_dat.txt","w") do io
    for i in 26:31
        print(io,L"{\footnotesize $f_{", i, L"}$}" )

        # Iterations
        xf_v = []
        q_v = []
        for k in 1:length(g_list)
            ds = setup_iterator(F_list[i], g_list[k], Float64.(X0[i]); algtype = :Steffensen)
            n, xf, q = compute_figure(ds, ε, max_it)
            push!(xf_v, xf)
            push!(q_v, q)
            # it = n.value
            print(io," & ", n)
        end
        
        #  Final point 
         for k in 1:length(g_list)
             if length(xf_v[k]) > 1
                 print(io," & (")
                 for x in xf_v[k]; print(io, round(Float64(x), digits =2), ", "); end
                 print(io,")")
             else
              print(io, " & ", round(Float64(xf_v[1]), digits =2)," ");
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



