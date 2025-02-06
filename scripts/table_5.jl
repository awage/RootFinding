using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function iterate(ds, x, ε, max_it)
    set_state!(ds, x) 
    n = @timed _get_iterations!(ds, ε, max_it)
    xf, _ = get_state(ds) 
    return n, xf
end

function print_table_all()
    ε = 1.e-14;  max_it = 50; force = false; Nsamples = Int(5e4)
    setprecision(BigFloat, 50; base = 10)

    open("table5_dat.txt","w") do io
    for i in  1:20
        print(io,L"{\footnotesize $f_{", i, L"}$}" )
        grid = ntuple(i -> range(-2, 2, length = 10), 1)
        for k in 1:length(g_list)
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

        ex_t = zeros(length(fam_list_real))
            N1 = stephenson_map_real(F_list[i], fam_list_real[1])
            ds1 = FunIterator(N1, [big(rand())])
            N2 = stephenson_map_real(F_list[i], fam_list_real[2])
            ds2 = FunIterator(N2, [big(rand())])
            N3 = stephenson_map_real(F_list[i], fam_list_real[3])
            ds3 = FunIterator(N3, [big(rand())])
            # ex_t[k] = exec_time
            t1 = 0.; t2 = 0.; t3 = 0.; k = 0
        while k < 100  
            x0 = big(4*(rand() - 0.5))
            n3, _ = iterate(ds3, [x0], ε, max_it)  
            if n3.value < max_it
                n2, _ = iterate(ds2, [x0], ε, max_it)  
                n1, _ = iterate(ds1, [x0], ε, max_it)  
                t1 += n1.time
                t2 += n2.time
                t3 += n3.time
                # @show n1.value, n2.value, n3.value
                # @show n1.time, n2.time, n3.time
                k = k + 1
            end
        end

        print(io," & ",  round(Float64(t1/t3), digits =2))
        print(io," & ",  round(Float64(t2/t3), digits =2))
        print(io," & ",  round(Float64(1.), digits =2))
        # for k in 1:length(fam_list_real)
        #     print(io," & ",  round(Float64(ex_t[k]/ex_t[3]), digits =2))
        # end
        
        println(io," \\\\")
    end

    ## Higher dimension functions


    for i in [1; 3:7]
        print(io,L"{\footnotesize $F_", i, L"$}")
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
        ex_t = zeros(length(fam_list_real))
        for k in 1:length(fam_list_real)
            N = stephenson_map_ndim(F2_list[i], fam_list_real[k], length(F2_X0[i]))
            d = _get_stats(N, F2_list[i], Nsamples, grid, ε, max_it; prefix = string("stats_F2_", i, "_g",k ), force = false)
            @unpack  exec_time = d
            ex_t[k] = exec_time
        end
        for k in 1:length(fam_list_real)
            print(io," & ",  round(Float64(ex_t[k]/ex_t[3]), digits =2))
        end
        println(io," \\\\")
        end
    end
end

print_table_all()



