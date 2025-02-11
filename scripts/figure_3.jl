using DrWatson
@quickactivate
using CairoMakie
using CodecZlib
using LaTeXStrings
using Statistics:mean
include("../src/function_stuff.jl")
include("../src/function_list.jl")
include("../src/color_stuff.jl")
include("../src/basins_compute.jl")

function plot_basins(N, res, ε = 1e-14, max_it =50; force = false, shaded = true, show_attractors = false, prefix = "stephenson")
    data = _get_basins(N, res, ε, max_it; force, prefix)
    @unpack basins, iterations, attractors, grid = data
    @show bas_num = unique(basins)
    if  1 < length(bas_num) < 100
        fig = plot_heatmap(grid, basins, iterations, attractors; ukeys = bas_num, shaded, show_attractors, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
        s = plotsdir(savename(string(prefix), @dict(res,ε),"png"))
        save(s, fig)
    end
end


# Plot all basins 
res = 300; 
for i in 1:20
for k in 1:length(g_list)
    F = [ x -> real(F_list[i](x[1]+im*x[2])),  x -> imag(F_list[i](x[1]+im*x[2]))]
    N = stephenson_map(F, g_list[k], 2)
    plot_basins(N, res; prefix = string("stephenson_f",i, "_g", k), force = true, shaded = true)
end
end


