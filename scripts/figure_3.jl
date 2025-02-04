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
    if  1 < length(bas_num) < 100
        fig = plot_heatmap(grid, basins, iterations, attractors; ukeys = bas_num, shaded, show_attractors, xticksvisible = false, yticksvisible = false, xticklabelsvisible = false, yticklabelsvisible = false)
        s = plotsdir(savename(string(prefix), @dict(res,ε),"png"))
        save(s, fig)
    end
end


# Plot all basins 
res = 500; 
for k in 1:length(fam_list_real)
    N = stephenson_map_ndim(F2_list[7], fam_list_real[k], length(F2_X0[7]))
        plot_basins(N, res; prefix = string("stephenson_f",7, "_g", k), force = true, shaded = true)
end


