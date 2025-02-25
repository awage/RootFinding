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
    return n.value[1], n.time
end

function print_table_all()
    ε = 1.e-8;  max_it = 100; force = true; Nsamples = Int(1e4)

    open("table5_dat.txt","w") do io
    for i in 1:20
        print(io,L"{\footnotesize $f_{", i, L"}$}" )
        grid = ntuple(i -> range(-1, 1, length = 10), length(X0[i]))
        it = zeros(length(g_list))
        ex = zeros(length(g_list))
        for k in 1:length(g_list)
            ds = setup_iterator(F_list[i], g_list[k], Float64.(X0[i]); algtype = :Steffensen)
            d = _get_stats(ds, Nsamples, grid, ε, max_it; prefix = string("stats_f", i, "_g",k ), force = force)
            @unpack nc, iterations, exec_time = d
            it[k] = iterations; ex[k] = exec_time
            @show nc, it[k], ex[k]
            print(io," & ",  round(Float64(nc*100), digits =1))
        end

        println(io," \\\\")

        for k in 1:length(g_list)
            print(io," & ",  round(Float64(it[k]), digits =1))
        end

        # for k in 1:length(g_list)
        #     print(io," & ",  round(Float64(ex[k]/ex[3]), digits =1))
        # end

        # N = [stephenson_map(F_list[i], g, length(X0[i]))  for g in g_list]
        # ds = [FunIterator(n, F_list[i], Float64.(X0[i])) for n in N] 
        # t1 = 0.; t2 = 0.; t3 = 0.; t4 = 0.; k = 0; cnt = 0;
        # sampler, = statespace_sampler(grid)
        # while k < 500  && cnt < Int(1e4)
        #     x0 = sampler()
        #     n4, dt4 = iterate(ds[4], x0, ε, max_it)  
        #     if n4 < max_it
        #         n3, dt3 = iterate(ds[3], x0, ε, max_it)  
        #         n2, dt2 = iterate(ds[2], x0, ε, max_it)  
        #         n1, dt1 = iterate(ds[1], x0, ε, max_it)  
        #         t1 += dt1
        #         t2 += dt2
        #         t3 += dt3
        #         t4 += dt4
        #         k = k + 1
        #     end
        #     cnt = cnt + 1
        # end

        # print(io," & ",  round(Float64(t1/t4), digits =2))
        # print(io," & ",  round(Float64(t2/t4), digits =2))
        # print(io," & ",  round(Float64(t3/t4), digits =2))
        # print(io," & ",  round(Float64(1.), digits =2))
        println(io," \\\\")
    end
end
end

print_table_all()



