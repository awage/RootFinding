using DrWatson
@quickactivate
using CodecZlib
using CairoMakie
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))



function print_table_all()
    res = 1000; ε = 1.e-14;  max_it = 50
    open("table2_dat.txt","w") do io
    for i in  [1:13; 15; 17]
        print(io, "{\\footnotesize ", string_list[i], "} " )
        N_β = beta_map(func_list[i])
        β4(x,y) =  2*x/(x+y)
        N_β_anneal = beta_map_anneal(func_list[i])


        # Mean iterations
        m0 = _get_mean_it(N_β, 0., i, res, ε, max_it)
        m1 = _get_mean_it(N_β, 1., i, res, ε, max_it)
        m_anneal = _get_mean_it(N_β_anneal, β4, i, res, ε, max_it; prefix = string("basins_anneal_", i))
        vec_col = set_color_cell([m0, m1, m_anneal])
        print(io, " & ", vec_col[1],  round(m0, digits =1))
        print(io, " & ", vec_col[2],  round(m1, digits =1))
        print(io, " & ", vec_col[3],  round(m_anneal, digits =1))

        # Non converging points 
        nc0 = _get_mean_nc(N_β, 0., i, res, ε, max_it)
        nc1 = _get_mean_nc(N_β, 1., i, res, ε, max_it)
        nc_anneal = _get_mean_nc(N_β_anneal, β4, i, res, ε, max_it; prefix = string("basins_anneal_", i))
        vec_col = set_color_cell([nc0, nc1, nc_anneal])
        print(io, " & ", vec_col[1],  round(Int, nc0*100))
        print(io, " & ", vec_col[2],  round(Int, nc1*100))
        print(io, " & ", vec_col[3],  round(Int, nc_anneal*100))

        # Computational time 
        t0_ref = _get_mean_t0(N_β, 0., i, res, ε, max_it)
        t1 = _get_mean_t0(N_β, 1., i, res, ε, max_it)
        t_anneal = _get_mean_t0(N_β_anneal, β4, i, res, ε, max_it; prefix = string("basins_anneal_", i))
        vec_col = set_color_cell([t0/t0_ref, t1/t0_ref, t_anneal/t0_ref])
        print(io, " & ", vec_col[1],  round(t0_ref/t0_ref, digits =2))
        print(io, " & ", vec_col[2],  round(t1/t0_ref, digits =2))
        print(io, " & ", vec_col[3],  round(t_anneal/t0_ref, digits =2))
        println(io," \\\\")
    end
    end
end

print_table_all()



