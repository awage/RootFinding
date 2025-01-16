using DrWatson
@quickactivate
using CairoMakie
using CodecZlib
using LaTeXStrings
using Statistics:mean
include("../src/function_stuff.jl")
include("../src/color_stuff.jl")
include("../src/basins_compute.jl")


function plot_basins(N, res, ε = 1e-14, max_it =50; force = false, shaded = true, show_attractors = false, prefix = "stephenson")
    data = _get_basins(N, res, ε, max_it; force, prefix)
    @unpack basins, iterations, attractors, grid = data
    @show bas_num = unique(basins)
    if length(bas_num) > 1
        fig = plot_heatmap(grid, basins, iterations, attractors; ukeys = bas_num, shaded, show_attractors, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
        s = plotsdir(savename(string(prefix), @dict(res,ε),"png"))
        save(s, fig)
    end
end


# Plot all basins 
res = 200; 
for (k,f) in enumerate(func_list)
    N = stephenson_map(func_list[k])
    plot_basins(N, res; prefix = string("stephenson_f",k), force = true)
end 
# res = 200
# f(z) = z^3 - z
# N = stephenson_map(f)
# data = _get_basins(N, 0, 0, res, 1e-8, 30; force = true)
# @unpack basins = data
# ds = DiscreteDynamicalSystem(N, [0.1, 0.2])
# z,t = trajectory(ds, 100, rand(2))

