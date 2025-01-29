using DrWatson
@quickactivate
using CodecZlib
# using CairoMakie
using LaTeXStrings
using Statistics
include(srcdir("function_stuff.jl"))
include(srcdir("basins_compute.jl"))

function compute_figure(N, x, ε, max_it)
    ds = DiscreteDynamicalSystem(N, [x])
        # set_state!(ds, [x])
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
    open("table3_dat.txt","w") do io
    for i in  1:5
            println(string_list_benchmark[i])
            print(io,"{\\footnotesize ", string_list_benchmark[i], "}" )

        # Iterations
        for k in 1:length(fam_list)
            N = stephenson_map_real(func_list_benchmark[i], fam_list_real[k])
            x0 = big(x0_list_benchmark[k])
            n, xf = compute_figure(N, x0, ε, max_it)
            it = n.value
            print(io," & ",  round(it, digits =1))
        end
        
        # Final point 
        for k in 1:length(fam_list)
            N = stephenson_map_real(func_list_benchmark[i], fam_list_real[k])
            x0 = big(x0_list_benchmark[k])
            n, xf = compute_figure(N, x0, ε, max_it)
            @show Float64(xf[1])
            # df = log10(abs(xf[1] - root_list_benchmark[k]))
            print(io," & ",  round(Float64(xf[1]), digits =1))
        end
        println(io," \\\\")
        end
    end
end

print_table_all()



