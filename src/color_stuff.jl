using Random: shuffle!, Xoshiro
using Colors
using Statistics:median

function colors_from_keys(ukeys)
    # Unfortunately, until `to_color` works with `Cycled`,
    # we need to explicitly add here some default colors...
    COLORS =[ 
        :black,
        :red,
        :green, 
        :blue,
        :purple,
        :yellow,]
    n = length(COLORS)
    v_col = []
    vals = zeros(2*length(ukeys))
    for k in eachindex(ukeys)
        push!(v_col, COLORS[mod(k-1,n)+1])
    end

    return Dict(k => v_col[i] for (i, k) in enumerate(ukeys))
end

function markers_from_keys(ukeys)
    MARKERS = [:circle, :dtriangle, :rect, :star5, :xcross, :diamond,
    :hexagon, :cross, :pentagon, :ltriangle, :rtriangle, :hline, :vline, :star4,]
    markers = Dict(k => MARKERS[mod1(i, length(MARKERS))] for (i, k) in enumerate(ukeys))
    return markers
end


function custom_colormap(ukeys)
    # Unfortunately, until `to_color` works with `Cycled`,
    # we need to explicitly add here some default colors...
    COLORS =[ 
        :black,
        :red,
        :green, 
        :blue,
        :purple,
        :yellow,]
    LIGHT_COLORS = [
        :gray95,
        :lightsalmon,
        :darkseagreen1,
        :azure,
        :lavenderblush,
        :lightyellow,]      
    DARK_COLORS = [
        :black,
        :red4,
        :darkgreen,
        :navyblue,
        :purple4,
        :gold4,]      
    n = length(COLORS)
    v_col = []
    vals = zeros(2*length(ukeys))
    for k in eachindex(ukeys)
        # push!(v_col, COLORS[mod(k-1,n)+1])
        # push!(v_col, :black)
        push!(v_col, LIGHT_COLORS[mod(k-1,n)+1])
        push!(v_col, DARK_COLORS[mod(k-1,n)+1])
        vals[2*k-1] = k-1
        vals[2*k] = k-1+0.9999999999
    end
    return cgrad(v_col, vals/maximum(vals))
end

function plot_heatmap(grid, basins, iterations, attractors; ukeys = unique(basins), shaded = true, show_attractors = true, kwargs...) 
    sort!(ukeys) # necessary because colormap is ordered
    ids = 1:length(ukeys)
    replace_dict = Dict(k => i for (i, k) in enumerate(ukeys))
    basins_to_plot = replace(basins.*1., replace_dict...)
    access = SVector(1,2)
    
    cmap = custom_colormap(ukeys)
    colors = colors_from_keys(ukeys)
    markers = markers_from_keys(ukeys)
    labels = Dict(ukeys .=> ukeys)
    add_legend = length(ukeys) < 7

    if shaded 
        maxit = median(iterations)*2
        it = findall(iterations .> maxit)
        iterations[it] .= maxit
        for i in ids 
            ind = findall(basins_to_plot .== i) 
            mn = minimum(iterations[ind])
            mx = maximum(iterations[ind])
            basins_to_plot[ind] .= basins_to_plot[ind] .+ 0.99.*(iterations[ind].-mn)/mx
        end
    end

    fig = Figure(size = (800,800))
    ax = Axis(fig[1,1]; kwargs...)
    heatmap!(ax, grid..., basins_to_plot;
        colormap = cmap,
        colorrange = (ids[1], ids[end]+1),
    )
    # Scatter attractors
    if show_attractors 
        for (i, k) ∈ enumerate(ukeys)
            k ≤ 0 && continue
            A = attractors[k]
            x, y = columns(A[:, access])
            scatter!(ax, x, y;
                color = colors[k], markersize = 20,
                marker = markers[k],
                strokewidth = 1.5, strokecolor = :white,
                label = "$(labels[k])"
            )
        end
        # Add legend using colors only
        add_legend && axislegend(ax)
    end

  return fig
end

