using Attractors
using LinearAlgebra:norm
using ProgressMeter

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


function _get_stats(N, f, Nsamples, grid, ε, max_it; prefix = "stats_", force = false)
    d = @dict(N, f, Nsamples, ε, max_it, grid) # parametros
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

# This is where the iterations are computed until 
# the stopping criterion is met
function _get_iterations!(ds, f, ε, max_it)
    xn_1 = get_state(ds) 
    step!(ds)
    fx = length(xn_1) > 1 ? map(h -> h(xn_1), f) : f(xn_1[1])
    xn = get_state(ds) 
    k = 1
    # stopping criterion is ∥x_n - x_{n-1}∥ + ∥f(x_{n-1})∥ ≤ ε
    while norm(xn - xn_1) + norm(fx) > ε  
        (k > max_it) && break 
        xn_1 = xn
        try 
            fx = length(xn_1) > 1 ? map(h -> h(xn_1), f) : f(xn_1[1])
            step!(ds)
        catch 
            @show xn_1, fx
            k = max_it + 1 
            break 
        end
        xn = get_state(ds) 
        k += 1
    end
    return k
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
    ds = DiscreteDynamicalSystem(N, [0.1, 0.2])
    xg = yg = range(-10, 10; length = 20001)
    grid = (xg, yg)
    # We set up a mapper so that we can identify roots automatically  
    mapper_beta = AttractorsViaRecurrences(ds, (xg, yg);
            sparse = true, consecutive_recurrences = 3000
    )
    xg = yg = range(-1, 1; length = res)
    grid = (xg, yg)

    basins = zeros(Int32,res,res); iterations = zeros(Int16,res,res)
    exec_time = zeros(res,res)

@showprogress for (i,x) in enumerate(xg), (j,y) in enumerate(yg) 
        set_state!(ds, [x,y])
        n = @timed _get_iterations!(ds,ε,max_it)
        if n.value > max_it
            # the alg. did not converge
            basins[i,j] = -1
        else
            # We identify the root with the mapper.
            basins[i,j] = mapper_beta([x,y])
        end
        iterations[i,j] = n.value
        exec_time[i,j] = n.time
    end

    Sb, Sbb = basin_entropy(basins) 
    _,_,fdim = basins_fractal_dimension(basins)
    attractors = extract_attractors(mapper_beta)
     
    x,y = choose_valid_ic!(ds, max_it, ε) 
    q = estimate_ACOC!(ds, 200,ε, x, y)
    
    return @strdict(grid, basins, iterations, exec_time, attractors, Sb, Sbb, fdim, q)
end


"""
Compute stats!
"""
function compute_stats(d)
    @unpack N, f,  Nsamples, ε, max_it, grid = d
    dim = length(grid)
    ds = DiscreteDynamicalSystem(N, big.(rand(dim)))
    # iterations = zeros(Int16, Nsamples)
    # exec_time = zeros(Nsamples)
    iterations = 0.0
    exec_time = 0.0
    nc = 0

    sampler, = statespace_sampler(grid)
    
    for k in 1:Nsamples
        set_state!(ds, big.(sampler()))
        n = @timed _get_iterations!(ds, f, ε, max_it)
        if n.value > max_it
            # the alg. did not converge
            nc += 1
        else
            iterations += n.value
            exec_time += n.time
        end
    end
    exec_time = exec_time/iterations
    iterations = iterations/(Nsamples - nc)
    nc = nc/Nsamples
    # @show nc, iterations, exec_time
    return @strdict(grid, Nsamples, iterations, exec_time, nc)
end

function choose_valid_ic!(ds, max_it, ε) 
 # make sure we pick an IC that converge to a root
 # with enough iterations (at least 8). 
     x = 0.; y = 0.; k = 0
     while true 
         x = 4*(rand()-0.5)
         y = 4*(rand()-0.5)
         set_state!(ds, [x,y])
         n = _get_iterations!(ds,ε,max_it)
         if (n < max_it) && (n ≥ 10)
            break
         end
         (k < 1000) || break
         k = k + 1
     end
     return x,y
 end


# Estimate order
function  estimate_ACOC!(ds, T, ε, x, y)
    yy,t = trajectory(ds, T, [x,y])
    qn_1 = 10000
    qn = qn_1 - 1
    k = 3
    while norm(yy[k+1] - yy[k]) > ε
        (k > T-2) && break 
        num = log(norm(yy[k+1] - yy[k])) - log(norm(yy[k] - yy[k-1])) 
        den = log(norm(yy[k] - yy[k-1]))- log(norm(yy[k-1] - yy[k-2]))
        qn_1 = qn
        qn = num/den
        k = k + 1 
    end
    return qn
end

function _get_mean_it(f, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    mit = mean(iterations[ind])
    return  mit
end

function _get_mean_t0(f, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    t0_ref = mean(exec_time[ind])
    return  t0_ref
end

function _get_mean_nc(f, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    nc = 1-length(ind)/length(basins)
    return  nc
end

function _get_mean_ps(f, i, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    ps_ref = length(ind)/sum(exec_time[ind])
    return  ps_ref
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
