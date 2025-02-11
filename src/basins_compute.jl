using Attractors
using LinearAlgebra:norm
using ProgressMeter

include(srcdir("function_stuff.jl"))


""" 
    function _get_basins(N,i,res,ε,max_it) -> data

Convenience function to compute and store the basins
and attractors of the funcion i. with the proximity algorithm

"""
function _get_basins(N, res, ε, max_it; prefix = "basins_", force = false)
    d = @dict(N, res, ε, max_it) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        compute_basins, # function
        prefix = prefix, # prefix for savename
        force = force, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    return data
end


function _get_stats(N, Nsamples, grid, ε, max_it; prefix = "stats_", force = false)
    d = @dict(N, Nsamples, ε, max_it, grid) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        compute_stats, # function
        prefix = prefix, # prefix for savename
        force = force, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    return data
end



"""
    compute_basins(d) -> Dict

Compute the basins using first AttractorsViaRecurrences to 
locate the attractors (roots) of the function N_β. These 
attractors are passed to AttractorsViaProximity with a 
convergence criterion ε such that the algorithm stops 
when |f(x) - r| < ε. 
The basins, the iteration matrix, the metrics and the attractors
are returned into a name dictionnary. 
"""
function compute_basins(d)
    @unpack N,  res, ε, max_it = d
    f = function(x,p,t); y,_ = N(x) ; return SVector{2}(y) end
    ds = DiscreteDynamicalSystem(f, big.(rand(2)))
    di = FunIterator(N, big.(rand(2)))
    xg = yg = range(-10, 10; length = 20001)
    grid = (xg, yg)
    # We set up a mapper so that we can identify roots automatically  
    mapper_beta = AttractorsViaRecurrences(ds, (xg, yg);
            sparse = true, consecutive_recurrences = 3000
    )
    xg = yg = range(-2, 2; length = res)
    grid = (xg, yg)

    basins = zeros(Int32,res,res); iterations = zeros(Int16,res,res)
    exec_time = zeros(res,res)

@showprogress for (i,x) in enumerate(xg), (j,y) in enumerate(yg) 
        set_state!(di, big.([x,y]))
        n = @timed _get_iterations!(di, ε, max_it)
        it = n.value[1]
        if it > max_it
            # the alg. did not converge
            basins[i,j] = -1
        else
            # We identify the root with the mapper.
            xf, _ = get_state(di)
            basins[i,j] = mapper_beta(xf)
        end
        iterations[i,j] = it
        exec_time[i,j] = n.time
    end

    Sb, Sbb = basin_entropy(basins) 
    _,_,fdim = basins_fractal_dimension(basins)
    attractors = extract_attractors(mapper_beta)
     
    # x,y = choose_valid_ic!(ds, max_it, ε) 
    # q = estimate_ACOC!(ds, 200,ε, x, y)
    q = 2    
    return @strdict(grid, basins, iterations, exec_time, attractors, Sb, Sbb, fdim, q)
end


"""
Compute stats!
"""
function compute_stats(d)
    @unpack N, Nsamples, ε, max_it, grid = d
    dim = length(grid)
    ds = FunIterator(N, big.(dim == 1 ? rand() : rand(dim)))
    iterations = 0.0
    exec_time = 0.0
    nc = 0

    sampler, = statespace_sampler(grid)
    if dim == 1
        samp = () -> big(sampler()[1])
    else
        samp = () -> big.(sampler())
    end

    for k in 1:Nsamples
        set_state!(ds, samp())
        n = @timed _get_iterations!(ds, ε, max_it)
        it = n.value[1]
        if it > max_it
            # the alg. did not converge
            nc += 1
        else
            iterations += it
            exec_time += n.time
        end
    end
    exec_time = exec_time/iterations
    iterations = iterations/(Nsamples - nc)
    nc = nc/Nsamples
    return @strdict(grid, Nsamples, iterations, exec_time, nc)
end

function choose_valid_ic!(ds, max_it, ε, sampler) 
 # make sure we pick an IC that converge to a root
 # with enough iterations (at least 8). 
     x = 0.;  k = 0
     while true 
         x = sampler()
         set_state!(ds, x)
         n = _get_iterations!(ds,ε,max_it)
         if (n < max_it) && (n ≥ 10)
            break
         end
         (k < 1000) || break
         k = k + 1
     end
     return x 
 end


# Estimate order
function  estimate_ACOC!(ds, T, yy)
    for k in 3:T-2
        num = log(norm(yy[k+1] - yy[k])) - log(norm(yy[k] - yy[k-1])) 
        den = log(norm(yy[k] - yy[k-1]))- log(norm(yy[k-1] - yy[k-2]))
        qn = num/den
    end
    return qn
end


function _get_q(N, res, ε, max_it; kwargs...)
    ds = DiscreteDynamicalSystem(N, [0.1, 0.2])
    x,y = choose_valid_ic!(ds, max_it, ε) 
    q = estimate_ACOC!(ds, 200,ε, x, y)
    return  q
end


# This small functions sets a color gradient between red and green depending 
# on the values in the input array
function set_color_cell(v) 
    mx = maximum(v)
    mn = minimum(v) 
    cl = @. round(Int, 100*(v/(mx-mn)) - (mn*(100/(mx-mn))))
    vc = [string("\\cellcolor{red!", c , "!green!15}") for c in cl]
    return vc 
end
