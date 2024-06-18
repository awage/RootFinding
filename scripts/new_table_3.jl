using DrWatson
@quickactivate
using CodecZlib
using CairoMakie
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))
include(srcdir("estimate_order.jl"))

function _get_mean_it(f, β, i, res, ε, max_it)
    data0 = _get_basins(f, β, i, res, ε, max_it)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    mit = mean(iterations[ind])
    return  mit
end

function _get_mean_t0(f, β, i, res, ε, max_it)
    data0 = _get_basins(f, β, i, res, ε, max_it)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    t0_ref = mean(exec_time[ind])
    return  t0_ref
end

function _get_mean_nc(f, β, i, res, ε, max_it)
    data0 = _get_basins(f, β, i, res, ε, max_it)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    nc = 1-length(ind)/length(basins)
    return  nc
end

function _get_mean_ps(f, β, i, res, ε, max_it)
    data0 = _get_basins(f, β, i, res, ε, max_it)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    ps_ref = length(ind)/sum(exec_time[ind])
    return  ps_ref
end

function _get_q(f, β, i, res, ε, max_it)
    data0 = _get_basins(f, β, i, res, ε, max_it)
    @unpack q = data0
    return  q
end

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
    # β_range = range(-1,1, step  = 0.5)
    res = 1000; ε = 1.e-15; Nsim = 1000; T = 1000; max_it = 50

    for i in  [1:13; 15; 17]
        print(string_list[i] )
        print(" ")
        N_β = beta_map(func_list[i])
        N_β_anneal = beta_map_anneal(func_list[j]); 

        t0_ref = _get_mean_t0(N_β, 0., i, res, ε, max_it)
        ps_ref = _get_mean_ps(N_β, 0., i, res, ε, max_it)

        # Mean iterations
        m0 = _get_mean_it(N_β, 0, i, res, ε, max_it)
        m1 = _get_mean_it(N_β, 1, i, res, ε, max_it)
        m_anneal = _get_mean_it(N_β_anneal, 1, i, res, ε, max_it; preffix = string("basins_prox_", i))
        vec_col = set_color_cell(m_it)
        for k in eachindex(β_range)
            print(" & ", vec_col[k],  round(m_it[k], digits =1))
        end

        # Non converging points 
        nc = [_get_mean_nc(N_β, β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(nc)
        for k in eachindex(β_range)
            print(" & ", vec_col[k],  round(Int,100*nc[k]))
        end

        # Computational time 
        t0 = [_get_mean_t0(N_β, β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(t0./t0_ref)
        for k in eachindex(β_range)
            print(" & ", vec_col[k],  round(t0[k]/t0_ref, digits =2))
        end
        println(" \\\\")
    end
        # println("\\end{tabular}")
        # println(" ")
end

print_table_all()



