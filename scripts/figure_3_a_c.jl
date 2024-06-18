using DrWatson
@quickactivate
using CairoMakie
using CodecZlib
using LaTeXStrings
using Statistics
include("../src/function_stuff.jl")
include("../src/basins_compute.jl")


function plot_f7_β()
# Plot metrics as a function of β
    β_range = range(0.0,1,step = 0.01)
    res = 300; ε = 1e-15; max_it = 50
    f = func_list[7]; i = 7
    Nit = zeros(length(β_range))
    stdNit = zeros(length(β_range))
    for (k,β) in enumerate(β_range)
        data = _get_basins(f, β, i, res, ε, max_it)
        @unpack iterations, basins, fdim = data
        ind = findall(basins .!= -1)
        Nit[k] = mean(iterations[ind])
        stdNit[k] = std(iterations[ind])
    end

    f = Figure(size = (600,500))
# gb = f[1,1] = GridLayout()
    ax = Axis(f[1,1], ylabel = L"N_{iter}", xlabel = L"\beta", yticklabelsize = 30, xticklabelsize = 30, ylabelsize = 40, xlabelsize = 40,  titlesize = 30 )
    band!(ax, β_range, Nit-stdNit, Nit+stdNit; color = (:lightblue,0.5))
    lines!(ax, β_range, Nit, color = :black)
    s = plotsdir(string("fig_Nit_f", i, ".png"))
    save(s, f)
end

function plot_f7_ε()
    ε_range = 10 .^range(-10,-2, length = 30)
    res = 300; max_it = 50
    m0 = zeros(length(ε_range))
    m1 = zeros(length(ε_range))
    f = func_list[7]; i = 7 
    for (k,ε) in enumerate(ε_range)
        data0 = _get_basins(f, 0, i, res, ε, max_it)
        data1 = _get_basins(f, 1, i, res, ε, max_it)
        @unpack iterations,basins = data0
        ind = findall(basins .!= -1)
        m0[k] = mean(iterations[ind])
        @unpack iterations,basins = data1
        ind = findall(basins .!= -1)
        m1[k] = mean(iterations[ind])
    end

    f = Figure(size = (600,500))
    ax = Axis(f[1,1],  ylabel = L"N_{iter}", xlabel = L"\varepsilon", yticklabelsize = 30, xticklabelsize = 30, ylabelsize = 40, xlabelsize = 40,  titlesize = 30, xscale = log10)
    lines!(ax, ε_range, m0, color = :black, label = L"\beta = 0")
    lines!(ax, ε_range, m1, color = :red, label = L"\beta = 1")
    s = plotsdir(string("fig_convergence_f", i, ".png"))
    axislegend(ax; labelsize = 30);
    save(plotsdir(s),f)
end


plot_f7_β()
plot_f7_ε()
