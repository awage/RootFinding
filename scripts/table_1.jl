using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Statistics
using Roots
using LinearAlgebra:norm
include("../src/function_stuff.jl")
include("../src/basins_compute.jl")

function _estimate_ACOC!(ds,Nsim)
    q_conv = zeros(Nsim)
    for n in 1:Nsim
        # pick random attractor: 
        rt = rand(roots)
        if !(i == 5 && rt ≈ -1.5)
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
        else 
            n -= 1 # back a step
        end
    end
    return q_conv
end

function estimate_q_real(i, β_vec, Nsim, T, ε) 
    q_conv = zeros(Nsim,length(β_vec))
    mean_q = zeros(length(β_vec))
    var_q = zeros(length(β_vec))
    roots = find_zeros(func_list[i], -10, 10; atol = ε)
    for (m,β) in enumerate(β_vec)
        @show β
        Nβ! = beta_map_real(func_list[i])
        ds = DiscreteDynamicalSystem(Nβ!, [0.2], [β])
        q_conv[:,m] = _estimate_ACOC!(ds,Nsim)
        ind = findall( 0 .< q_conv[:,m] .< 10)
        mean_q[m] = mean(q_conv[ind,m])
        var_q[m] = var(q_conv[ind,m])
    end
    return mean_q,var_q, q_conv
end

function estimate_q_cmplx(i, β_vec, Nsim, T, ε = 1e-5)
    q_conv = zeros(Nsim,length(β_vec))
    mean_q = zeros(length(β_vec))
    var_q = zeros(length(β_vec))
    for (m,β) in enumerate(β_vec)
        @show β
        data0 = _get_dat(func_list[i], β, i, 100, ε)
        @unpack basins, attractors, grid = data0
        Nβ = beta_map(func_list[i])
        ds = DiscreteDynamicalSystem(Nβ, [0.1, 0.2], [β])
        q_conv[:,m] = _estimate_ACOC!(ds,Nsim)
        ind = findall( 0 .< q_conv[:,m] .< 10)
        mean_q[m] = mean(q_conv[ind,m])
        var_q[m] = var(q_conv[ind,m])
    end
    return mean_q,var_q, q_conv
end

# Estimation of convergence order in the complex plane: 
function get_order_functions_cmplx()
    Nsim = 5000
    T = 1000
    β_vec = [0., 1.]
    m = zeros(15,length(β_vec))
    v = zeros(15,length(β_vec))
    for i in [1:13; 15]
        m[i,:],v[i,:],_ = estimate_q_cmplx(i, β_vec, Nsim, T, 1e-15) 
    end

    for k in [1:13; 15]
        println("function ", k, " | ", round(m[k,1], digits =2), " ±", round(v[k,1],digits =2), " | ", round(m[k,2], digits =2), " ±", round(sqrt(v[k,2]), digits =2))
    end
end


# Estimation of convergence order on the real line: 
function get_order_functions_real()
    Nsim = 5000
    T = 1000
    # β_vec = range(0,1, step= 0.1)
    β_vec = [0., 1.]
    m = zeros(16,length(β_vec))
    v = zeros(16,length(β_vec))
    for i in 1:14
        m[i,:],v[i,:],_ = estimate_q_real(i, β_vec, Nsim, T, 1e-15) 
    end

    for i in 1:14 
        println(string_list[i], " & ", round(m[i,1], digits =2),   " & ", round(m[i,2], digits =2) , "\\\\" )
    end
    return m,v
end


# get_order_functions_cmplx()
get_order_functions_real()
