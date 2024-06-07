using Attractors
using LinearAlgebra:norm

""" 
    function _get_basins(f,β,i,res,ε) -> data

Convenience function to compute and store the basins
and attractors of the funcion i. with the proximity algorithm

"""
function _get_basins(f,β,i,res,ε)
    N_β = beta_map(f)
    d = @dict(N_β, β, res, ε) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        compute_basins_prox, # function
        prefix = string("basins_prox_",i), # prefix for savename
        force = false, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    return data
end

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
    compute_basins_prox(d) -> Dict

Compute the basins using first AttractorsViaRecurrences to 
locate the attractors (roots) of the function N_β. These 
attractors are passed to AttractorsViaProximity with a 
convergence criterion ε such that the algorithm stops 
when |f(x) - r| < ε. 
The basins, the iteration matrix, the metrics and the attractors
are returned into a name dictionnary. 
"""
function compute_basins_prox(d)
    @unpack N_β, β, res, ε, max_it = d
    ds = DiscreteDynamicalSystem(N_β, [0.1, 0.2], [β])
    # Use non-sparse for using `basins_of_attraction`
    xg = yg = range(-10, 10; length = 10000)
    grid = (xg, yg)
    mapper_beta = AttractorsViaRecurrences(ds, (xg, yg);
            sparse = true, consecutive_recurrences = 3000
    )
    xg = yg = range(-2, 2; length = res)
    grid = (xg, yg)

    basins = zeros(Int8, res,res); iterations = zeros(Int16,res,res)
    exec_time = zeros(res,res)
    # if there is not enough attractors return empty matrix
    # if length(attractors) < 2
    #     Sb = Sbb = fdim = 0.
    #     return @strdict(β, grid, basins, iterations, attractors, Sb, Sbb, fdim)
    # end

    # Proceed to compute the basins and the iterations.
    # mapper_prox = AttractorsViaProximity(ds, attractors, ε; Ttr = 0, horizon_limit = 1e50, consecutive_lost_steps = 10000) 
    
    for (i,x) in enumerate(xg), (j,y) in enumerate(yg) 
        set_state!(ds, [x,y])
        n = @timed _get_iterations!(ds,ε,max_it)
        if n.value > max_it
            basins[i,j] = -1
        else
            basins[i,j] = mapper_beta([x,y])
        end
        iterations[i,j] = n.value
        exec_time[i,j] = n.time
    end

    Sb, Sbb = basin_entropy(basins) 
    _,_,fdim = basins_fractal_dimension(basins)
    return @strdict(β, grid, basins, iterations, exec_time, attractors, Sb, Sbb, fdim)
end



