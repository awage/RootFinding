using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Statistics
include("../src/function_stuff.jl")
include("../src/basins_compute.jl")
include(srcdir("estimate_order.jl"))


function print_table_all()
    β_range = range(-1,1, step  = 0.5)
    res = 300; ε = 1.e-8; Nsim = 1000; T = 1000
    max_it = 50
    mit= zeros(15,length(β_range))
    t0 = zeros(15,length(β_range))
    nc = zeros(15,length(β_range))
    q = zeros(15,length(β_range))
    ps = zeros(15,length(β_range))
    its = zeros(15,length(β_range))
    ff = Figure(size = (600,500))
    ax = Axis(ff[1,1], 
    ylabel = L"\eta", xlabel = L"\beta", yticklabelsize = 30, xticklabelsize = 30, ylabelsize = 40, xlabelsize = 40,  titlesize = 30)
    for i in [1:13; 15]
        println("|funct: ", i, " |  IT  | EX μs| NC %  | OR |")
        for (k,β) in enumerate(β_range)
            data0 = _get_basins(func_list[i], β, i, res, ε)
            @unpack attractors,iterations,basins,exec_time = data0
            roots = [attractors[kk][1] for kk in 1:length(attractors)]
            q, qv = get_ACOC(roots, 5000, i, β, 1000, 1e-15)
            ind = findall(basins .!= -1)
            mit = mean(iterations[ind])
            t0 = mean(exec_time[ind])
            nc = 1-length(ind)/length(basins)
            ps = length(ind)/sum(exec_time[ind])
            its = sum(iterations[ind])/sum(exec_time[ind])
            
            println("| ", "β = ", β, 
                " | ", round(q, digits =2),
                " | ", round(mit, digits =2), 
                " | ", round(nc*100, digits =2),
                " | ", round(t0*1e6, digits =2), 
                " | ", round(ps*1e-6, digits =2),
                " | ", round(its*1e-6, digits =2))

        end
    end
end
print_table_all()
