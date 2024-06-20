using DrWatson
@quickactivate
using CodecZlib
using CairoMakie
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))
# include(srcdir("estimate_order.jl"))


function set_color_cell(v) 
    mx = maximum(v)
    mn = minimum(v) 
# int(round(100*(#1/(\maxval-\minval))-(\minval*(100/(\maxval-\minval)    
    cl = @. round(Int, 100*(v/(mx-mn)) - (mn*(100/(mx-mn))))
    # cl = (v .- mn)./mx
    vc = [string("\\cellcolor{red!", c , "!green!15}") for c in cl]
    return vc 
end

function print_table_all()
    β_range = range(-1,1, step  = 0.5)
    res = 1000; ε = 1.e-14;  max_it = 50

    open("table1_dat.txt","w") do io
    for i in  [1:13; 15; 17]
        println(string_list[i])
        print(io,"{\\footnotesize ", string_list[i], "}" )
        N_β = beta_map(func_list[i])

        t0_ref = _get_mean_t0(N_β, 0., i, res, ε, max_it)
        ps_ref = _get_mean_ps(N_β, 0., i, res, ε, max_it)

        # Mean iterations
        m_it = [_get_mean_it(N_β, β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(m_it)
        for k in eachindex(β_range)
            print(io," & ", vec_col[k],  round(m_it[k], digits =1))
        end

        # Non converging points 
        nc = [_get_mean_nc(N_β, β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(nc)
        for k in eachindex(β_range)
            print(io, " & ", vec_col[k],  round(Int,100*nc[k]))
        end

        # Computational time 
        t0 = [_get_mean_t0(N_β, β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(t0./t0_ref)
        for k in eachindex(β_range)
            print(io, " & ", vec_col[k],  round(t0[k]/t0_ref, digits =2))
        end
        println(io," \\\\")
    end
    end
end

print_table_all()



