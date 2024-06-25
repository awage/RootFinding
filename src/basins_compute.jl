using Attractors
using LinearAlgebra:norm
using ProgressMeter

""" 
    function _get_basins(N_β,β,i,res,ε,max_it) -> data

Convenience function to compute and store the basins
and attractors of the funcion i. with the proximity algorithm

"""
function _get_basins(N_β,β,i,res,ε,max_it; prefix = string("basins_prox_", i))
    # N_β = beta_map(f)
    d = @dict(N_β, β, res, ε, max_it) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        compute_basins, # function
        prefix = prefix, # prefix for savename
        force = false, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    return data
end

# This is where the iterations are computed until 
# the stopping criterion is met
function _get_iterations!(ds, ε, max_it)
    xn_1 = get_state(ds) 
    step!(ds)
    xn = get_state(ds) 
    k = 1
    # stopping criterion is ∥x_n - x_{n-1}∥ ≤ ε
    while norm(xn - xn_1) > ε
        (k > max_it) && break 
        xn_1 = xn
        step!(ds)
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
    @unpack N_β, β, res, ε, max_it = d
    ds = DiscreteDynamicalSystem(N_β, [0.1, 0.2], [β])
    xg = yg = range(-10, 10; length = 10000)
    grid = (xg, yg)
    # We set up a mapper so that we can identify roots automatically  
    mapper_beta = AttractorsViaRecurrences(ds, (xg, yg);
            sparse = true, consecutive_recurrences = 3000
    )
    xg = yg = range(-2, 2; length = res)
    grid = (xg, yg)

    basins = zeros(Int8, res,res); iterations = zeros(Int16,res,res)
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
    @show q = estimate_ACOC!(ds, 200,ε, x, y)
    

    return @strdict(β, grid, basins, iterations, exec_time, attractors, Sb, Sbb, fdim, q)
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

function _get_mean_it(f, β, i, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, β, i, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    mit = mean(iterations[ind])
    return  mit
end

function _get_mean_t0(f, β, i, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, β, i, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    t0_ref = mean(exec_time[ind])
    return  t0_ref
end

function _get_mean_nc(f, β, i, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, β, i, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    nc = 1-length(ind)/length(basins)
    return  nc
end

function _get_mean_ps(f, β, i, res, ε, max_it; kwargs...)
    data0 = _get_basins(f, β, i, res, ε, max_it; kwargs...)
    @unpack iterations,basins,  exec_time = data0
    ind = findall(basins .!= -1)
    ps_ref = length(ind)/sum(exec_time[ind])
    return  ps_ref
end

function _get_q(N_β, β, i, res, ε, max_it; kwargs...)
    ds = DiscreteDynamicalSystem(N_β, [0.1, 0.2], [β])
    x,y = choose_valid_ic!(ds, max_it, ε) 
    @show q = estimate_ACOC!(ds, 200,ε, x, y)
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
