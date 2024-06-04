
function _estimate_ACOC!(ds,Nsim)
    q_conv = zeros(Nsim)
    for n in 1:Nsim
        # pick random roots: 
        rt = rand(roots)
        yy,t = trajectory(ds, T, [rt .+ randn()*0.2])
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
            q_conv[n,m] = q[end]
        end
    end
    return q_conv
end

