using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function compute_figure(N, x, ε, max_it)
    ds = DiscreteDynamicalSystem(N, x)
        # set_state!(ds, [x])
        n = @timed _get_iterations!(ds,ε,max_it)
        xf = current_state(ds) 
        return n, xf
end

function print_table_all()
    ε = 1.e-14;  max_it = 50; force = true; Nsamples = Int(1e4)
    setprecision(BigFloat, 50; base = 10)

    open("table5_dat.txt","w") do io
    for i in  1:21
        print(io,L"{\footnotesize $f_{", i, L"}$}" )
        grid = ntuple(i -> range(-2, 2, length = 10), 1)
        for k in 1:length(fam_list)
            N = stephenson_map_real(F_list[i], fam_list_real[k])
            d = _get_stats(N, F_list[i], Nsamples, grid, ε, max_it; prefix = string("stats_f", i, "_g",k ), force = force)
            @unpack nc, iterations, exec_time = d
            @show nc, iterations, exec_time
            print(io," & ",  round(Float64(nc*100), digits =1))
        end

        for k in 1:length(fam_list)
            N = stephenson_map_real(F_list[i], fam_list_real[k])
            d = _get_stats(N, F_list[i], Nsamples,  grid, ε, max_it; prefix = string("stats_f", i, "_g",k ), force = false)
            @unpack nc, iterations, exec_time = d
            print(io," & ",  round(Float64(iterations), digits =1))
        end

        for k in 1:length(fam_list)
            N = stephenson_map_real(F_list[i], fam_list_real[k])
            d = _get_stats(N, F_list[i], Nsamples, grid, ε, max_it; prefix = string("stats_f", i, "_g",k ), force = false)
            @unpack nc, iterations, exec_time = d
            print(io," & ",  round(Float64(1e6*exec_time), digits =2))
        end
        println(io," \\\\")
    end

    ## Higher dimension functions


    for i in [1; 3:7]
        print(io,L"{\footnotesize $F_{", i, L"$}")
        grid = ntuple(i -> range(-2, 2, length = 10), length(F2_X0[i]))
        for k in 1:length(fam_list)
            N = stephenson_map_ndim(F2_list[i], fam_list_real[k], length(F2_X0[i]))
            d = _get_stats(N, F2_list[i], Nsamples, grid, ε, max_it; prefix = string("stats_F2_", i, "_g",k ), force = force)
            @unpack nc = d
            print(io," & ",  round(Float64(100*nc), digits =1))
        end
        for k in 1:length(fam_list)
            N = stephenson_map_ndim(F2_list[i], fam_list_real[k], length(F2_X0[i]))
            d = _get_stats(N, F2_list[i], Nsamples, grid, ε, max_it; prefix = string("stats_F2_", i, "_g",k ), force = false)
            @unpack  iterations = d
            print(io," & ",  round(Float64(iterations), digits =1))
        end
        for k in 1:length(fam_list)
            N = stephenson_map_ndim(F2_list[i], fam_list_real[k], length(F2_X0[i]))
            d = _get_stats(N, F2_list[i], Nsamples, grid, ε, max_it; prefix = string("stats_F2_", i, "_g",k ), force = false)
            @unpack  exec_time = d
            print(io," & ",  round(Float64(1e4*exec_time), digits =1))
        end
        println(io," \\\\")
        end
    end
end

print_table_all()



