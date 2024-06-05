function get_ACOC(roots, Nsim, i, β, T, ε)
    Nβ = beta_map(func_list[i])
    d = @dict(Nβ, β, Nsim, roots, ε, T) # parametros
    data, file = produce_or_load(
        datadir(""), # path
        d, # container for parameter
        _estimate_ACOC, # function
        prefix = string("convergence_",i), # prefix for savename
        force = false, # true for forcing sims
        wsave_kwargs = (;compress = true)
    )
    @unpack q_conv = data
    ind = findall( 0 .< q_conv .< 10)
    return mean(q_conv), var(q_conv)
end

function _estimate_ACOC(d)
    @unpack Nβ, β, Nsim, roots, ε, T = d
    ds = DiscreteDynamicalSystem(Nβ, [0.1, 0.2], [β])
    q_conv = zeros(Nsim)
    for n in 1:Nsim
        # pick random roots: 
        rt = rand(roots)
        ic =  rt .+ randn(length(roots[1]))*0.2
        yy,t = trajectory(ds, T, ic)
        v = [norm(yy[k] .- rt) for k in 1:T]
        indf = findfirst(v .< ε) 
        if !isnothing(indf) && indf > 3
            q = Vector{Float64}()
            for k in 3:(indf-1)
                num = log(norm(yy[k+1] - yy[k])/norm(yy[k] - yy[k-1])) 
                den = log(norm(yy[k] - yy[k-1])/norm(yy[k-1] - yy[k-2]))
                qe = num/den
                push!(q,qe)
            end
            q_conv[n] = q[end]
        end
    end
    return @strdict(q_conv)
end

