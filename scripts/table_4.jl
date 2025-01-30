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
    n = @timed _get_iterations!(ds,ε,max_it)
    xf = current_state(ds) 
    return n, xf
end

# Table preamble 
# \begin{tabular}{p{6cm} ddd | ddd |ddd}
 # & \multicolumn{3}{c}{Mean it. per pt} & \multicolumn{3}{c}{Non conv. (\%)} &\multicolumn{3}{c}{T per pt } \\
 # & \multicolumn{3}{c}{$\beta$} & \multicolumn{3}{c}{$\beta$} &\multicolumn{3}{c}{$\beta$} \\    
function print_table_all()
    ε = 1.e-14;  max_it = 50; 
    setprecision(BigFloat, 50; base = 10)

    open("table4_dat.txt","w") do io
    for i in  1:21
        print(io,"{\\footnotesize f}" )

        # Iterations
        x0 = big(F_X0[i])
        for k in 1:length(fam_list)
            N = stephenson_map_real(F_list[i], fam_list_real[k])
            n, xf = compute_figure(N, x0, ε, max_it)
            it = n.value
            print(io," & ",  round(it, digits =1))
        end
        
        # Final point 
        x0 = big(F_X0[i])
        for k in 1:length(fam_list)
            N = stephenson_map_real(F_list[i], fam_list_real[k])
            n, xf = compute_figure(N, x0, ε, max_it)
            @show Float64(xf[1])
            print(io," & ",  round(Float64(xf[1]), digits =1))
        end

        println(io," \\\\")
    end

    ## Higher dimension functions
    for i in  1:6
        # println(string_list_benchmark[i])
        print(io,"{\\footnotesize f}" )

        # Iterations
        @show x0 = BigFloat.(F2_X0[i])
        for k in 1:length(fam_list)
            N = stephenson_map_ndim(F2_list[i], fam_list_real[k], length(F2_X0[i]))
            n, xf = compute_figure(N, x0, ε, max_it)
            it = n.value
            print(io," & ",  round(it, digits =1))
        end
        
        # Final point 
        @show x0 = BigFloat.(F2_X0[i])
        for k in 1:length(fam_list)
        N = stephenson_map_ndim(F2_list[i], fam_list_real[k], length(F2_X0[i]))
            n, xf = compute_figure(N, x0, ε, max_it)
            @show Float64(xf[1])
            print(io," & ",  round.(Float64.(xf), digits =2))
        end

        println(io," \\\\")
        end
    end
end

print_table_all()



