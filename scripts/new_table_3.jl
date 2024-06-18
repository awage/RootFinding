using DrWatson
@quickactivate
using CodecZlib
using CairoMakie
using LaTeXStrings
using Statistics
include("../src/function_stuff.jl")
include("../src/basins_compute.jl")
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
    β_range = range(-1,1, step  = 0.5)
    # β_range = [0., -1, -0.5, 0.5, 1.]
    res = 1000; ε = 1.e-15; Nsim = 1000; T = 1000
    max_it = 50
    # mit= zeros(15,length(β_range))
    # t0 = zeros(15,length(β_range))
    # nc = zeros(15,length(β_range))
    # q = zeros(15,length(β_range))
    # ps = zeros(15,length(β_range))
    # its = zeros(15,length(β_range))
    ff = Figure(size = (600,500))
    ax = Axis(ff[1,1], 
    ylabel = L"\eta", xlabel = L"\beta", yticklabelsize = 30, xticklabelsize = 30, ylabelsize = 40, xlabelsize = 40,  titlesize = 30)
    for i in  [1:13; 15; 17]
        print(string_list[i] )
        print(" ")

        t0_ref = _get_mean_t0(func_list[i], 0., i, res, ε, max_it)
        ps_ref = _get_mean_ps(func_list[i], 0., i, res, ε, max_it)

        m_it = [_get_mean_it(func_list[i], β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(m_it)
        for k in eachindex(β_range)
            print(" & ", vec_col[k],  round(m_it[k], digits =1))
        end

        nc = [_get_mean_nc(func_list[i], β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(nc)
        for k in eachindex(β_range)
            print(" & ", vec_col[k],  round(Int,100*nc[k]))
        end

        t0 = [_get_mean_t0(func_list[i], β, i, res, ε, max_it) for β in β_range]
        vec_col = set_color_cell(t0./t0_ref)
        for k in eachindex(β_range)
            print(" & ", vec_col[k],  round(t0[k]/t0_ref, digits =2))
        end
        println(" \\\\")
    end
        println("\\end{tabular}")
        println(" ")
end
print_table_all()



