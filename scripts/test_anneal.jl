using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function compute_figure(ds, ε, max_it)
    n, yy = _get_iterations!(ds, ε, max_it)
    xf, fx = get_state(ds) 
    @show xf, fx
    if 5 ≤ n < max_it
        q = estimate_ACOC!(n, yy)
    else
        q = 0
    end
    return n, xf, q
end

ε = 1e-25;  max_it = 100; 
# ε = 1e-8;  max_it = 100; 
setprecision(BigFloat, 50; base = 10)


ds = setup_iterator(F_list[19], g_list[1], X0[19]; algtype = :Steffensen)
@show n, xf, q = compute_figure(ds, ε, max_it)

ds = setup_iterator(F_list[1], g_list[1], X0[1]; algtype = :Steffensen)
@show n, xf, q = compute_figure(ds, ε, max_it)


ds = setup_iterator(F_list[26], g_list[1], X0[26]; algtype = :Steffensen)
@show n, xf, q = compute_figure(ds, ε, max_it)

ds = setup_iterator(F_list[19], g_list[1], X0[19]; algtype = :accelerated)
@show n, xf, q = compute_figure(ds, ε, max_it)


ds = setup_iterator(F_list[1], g_list[1], X0[1]; algtype = :accelerated)
@show n, xf, q = compute_figure(ds, ε, max_it)

ds = setup_iterator(F_list[1], g_list[1], X0[1]; algtype = :accelerated_secant)
@show n, xf, q = compute_figure(ds, ε, max_it)

# ds = setup_iterator(F_list[26], g_list[1], X0[26]; algtype = :accelerated)
# @show n, xf, q = compute_figure(ds, ε, max_it)
