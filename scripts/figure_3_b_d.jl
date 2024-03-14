using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Statistics:mean
include("../src/function_stuff.jl")
include("../src/basins_compute.jl")


function plot_benchmark_β_all()
    β_range = range(0,1, step  = 0.05)
    res = 300; ε = 1.e-8 
    m0 = zeros(15,length(β_range))
        ff = Figure(size = (600,500))
        ax = Axis(ff[1,1], 
           ylabel = L"1 - \langle N_{\beta = 1}\rangle /\langle N_{\beta = 0}\rangle", xlabel = L"\beta", yticklabelsize = 30, xticklabelsize = 30, ylabelsize = 40, xlabelsize = 40,  titlesize = 30)

    for i in [1:13; 15]
        f = func_list[i]; 
        for (k,β) in enumerate(β_range)
            data0 = _get_dat(f, β, i, res, ε)
            @unpack iterations,basins = data0
            ind = findall(basins .!= -1)
            m0[i,k] = mean(iterations[ind])
        end
        @show mit = mean(m0[i,:]./m0[i,1])
        clr = mit > 1 ? :lightsalmon : :lightgreen
        lines!(ax, β_range, 1 .- m0[i,:]./m0[i,1], color = clr)
    end
        s = plotsdir("fig_bench_com_β_all.pdf")
        ylims!(ax, -1.5,0.5)
        save(plotsdir(s),ff)
end


function plot_benchmark_ε_all()
    res = 300; 
    ff = Figure(size = (600,500))
    ax = Axis(ff[1,1], 
    ylabel = L"1 - \langle N_{\beta = 1}\rangle /\langle N_{\beta = 0}\rangle", xlabel = L"f_k", yticklabelsize = 30, xticklabelsize = 20, ylabelsize = 40, xlabelsize = 40,  titlesize = 30, xticks = 1:15 )
    mv1 = zeros(15)
    mv0 = zeros(15)
    ev0 = zeros(15)
    ev1 = zeros(15)
    ε = 1e-9
    for i in [1:13 ; 15]
        f = func_list[i]; 
        data0 = _get_dat(f, 0, i, res, ε)
        data1 = _get_dat(f, 1, i, res, ε)
        @unpack iterations,basins = data0
        ind = findall(basins .!= -1)
        m0 = mean(iterations[ind])
        v0 = std(iterations[ind])
        @unpack iterations,basins = data1
        ind = findall(basins .!= -1)
        m1 = mean(iterations[ind])
        v1 = std(iterations[ind])
        mv1[i] = m1
        ev1[i] = v1
        mv0[i] = m0
        ev0[i] = v0
    end
        ind = 1:15
        # The values for i = 14 have been estimated 
        # with real values only
        mv1[14] = 24.01005586592179 
        mv0[14] = 40.092424242424244 
        barplot!(ax, ind, 1 .- (mv1[ind]./mv0[ind]); color = :lightsalmon)
        s = plotsdir("fig_benchmark_ε_all.png")
        save(plotsdir(s),ff)
end


plot_benchmark_ε_all()
plot_benchmark_β_all()
