using DrWatson
@quickactivate
using CodecZlib
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("function_list.jl"))
include(srcdir("basins_compute.jl"))

function iterate(ds, x, ε, max_it)
    d = length(x) 
    set_state!(ds, d == 1 ? x[1] : x) 
    n = @timed _get_iterations!(ds, ε, max_it)
    xf, _ = get_state(ds) 
    return n, xf
end

function print_table_all()
    ε = 1.e-14;  max_it = 50; force = false; Nsamples = Int(1e4)
    setprecision(BigFloat, 50; base = 10)

    open("table5_dat.txt","w") do io
    for i in  1:28
        print(io,L"{\footnotesize $f_{", i, L"}$}" )
        grid = ntuple(i -> range(-2, 2, length = 10), length(X0[i]))
        it = zeros(length(g_list))
        ex = zeros(length(g_list))
        for k in 1:length(g_list)
            N = stephenson_map(F_list[i], g_list[k], length(X0[i]))
            d = _get_stats(N, Nsamples, grid, ε, max_it; prefix = string("stats_f", i, "_g",k ), force = force)
            @unpack nc, iterations, exec_time = d
            it[k] = iterations; ex[k] = exec_time
            @show nc, it[k], ex[k]
            print(io," & ",  round(Float64(nc*100), digits =1))
        end

        for k in 1:length(g_list)
            print(io," & ",  round(Float64(it[k]), digits =1))
        end

        # for k in 1:length(g_list)
        #     print(io," & ",  round(Float64(ex[k]/ex[3]), digits =1))
        # end

        N = [stephenson_map(F_list[i], g, length(X0[i]))  for g in g_list]
        ds = [FunIterator(n, X0[i]) for n in N] 
        t1 = 0.; t2 = 0.; t3 = 0.; k = 0; cnt = 0;
        sampler, = statespace_sampler(grid)
        while k < 500  && cnt < Int(1e4)
            x0 = big.(sampler())
            n3, _ = iterate(ds[3], x0, ε, max_it)  
            if n3.value < max_it
                n2, _ = iterate(ds[2], x0, ε, max_it)  
                n1, _ = iterate(ds[1], x0, ε, max_it)  
                t1 += n1.time
                t2 += n2.time
                t3 += n3.time
                # @show n1.value, n2.value, n3.value
                # @show n1.time, n2.time, n3.time
                k = k + 1
            end
            cnt = cnt + 1
        end

        print(io," & ",  round(Float64(t1/t3), digits =2))
        print(io," & ",  round(Float64(t2/t3), digits =2))
        print(io," & ",  round(Float64(1.), digits =2))
        
        println(io," \\\\")
    end
end
end

print_table_all()



