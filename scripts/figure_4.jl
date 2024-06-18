using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Statistics:mean
include("../src/function_stuff.jl")
include("../src/color_stuff.jl")
include("../src/basins_compute.jl")


function plot_basins(f,β,i,res, ε = 1e-15, max_it =50; shaded = true, show_attractors = false)
    data = _get_basins(f, β, i, res, ε,max_it)
    @unpack basins, iterations, attractors, grid = data
    @show bas_num = unique(basins)

    if length(bas_num) > 1
        s = LaTeXString(string(string_list[i], " β = ", β))
        fig = plot_heatmap(grid, basins, iterations, attractors; ukeys = bas_num, shaded, show_attractors, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
        s = plotsdir(savename(string("fig_f", i), @dict(res,β,ε),"png"))
        save(s, fig)
    end
end


# Plot all basins 
β_range = [0.,1.]
f_list = [2,7,15]
res = 2000
for i in f_list, β in β_range
    plot_basins(func_list[i], β, i, res; shaded = true)
end

# Plot metrics as a function of β
β_range = range(-1.,1,step = 0.1)
res = 500; ε = 1e-15; max_it = 50
for i in f_list
    Sb_v = zeros(length(β_range))
    Sbb_v = zeros(length(β_range))
    fdim_v = zeros(length(β_range))
    for (k,β) in enumerate(β_range)
        data = _get_basins(func_list[i], β, i, res, ε,max_it)
        @unpack iterations, Sb, Sbb, fdim = data
        Sb_v[k],Sbb_v[k],fdim_v[k] = Sb, Sbb, fdim
    end
    f = Figure()
    ax = Axis(f[1,1], ylabel = L"S_b", xlabel = L"\beta")
    lines!(ax, β_range, Sb_v, color = :black)
    s = plotsdir(string("fig_metrics_prox_f", i, ".png"))
    save(s, f)
end

