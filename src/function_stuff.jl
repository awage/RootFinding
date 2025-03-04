using ForwardDiff: derivative
using LinearAlgebra

"""
Custum iterator, it returns the state but also the function 
evaluation. accepts also a variable beta for the accelerated 
method.
"""
mutable struct State{N <: Number}
    x::Union{Vector{N}, N}
    fx::Union{Vector{N}, N}
    dfx::Union{Matrix{N}, Vector{N}, N}
end

mutable struct FunIterator
    N!::Function
    f::Union{Function, Vector{Function}}
    S::State
end

function setup_iterator(f::Union{Function,Vector{Function}}, g::Function, x; algtype = :Steffensen)
    N! = if algtype == :Steffensen
        stephenson_map(f, g, length(x))
    elseif algtype == :accelerated
        stephenson_map_accel(f, g, length(x))
    elseif algtype == :accelerated_secant
        stephenson_map_accel_secant(f, g, length(x))
    else
        error("Invalid algtype: $algtype. Choose :Steffensen or :accelerated.")
    end
    
    fi = FunIterator(N!, f, State(zero(x),zero(x), zero(x)))  
    set_state!(fi, x)  # Initialize the state
    return fi
end


"""
Constructor for FunIterator.  Initializes the State struct.
"""
function FunIterator(N!::Function, f, x::Union{Vector{T}, T} where T <: Number; n = length(x), k = 1) 
    fx = isa(f, Vector) ? map(h -> h(x), f) : f(x)
    dfx = if isa(f, Vector)
        ones(eltype(x), n, length(f))
    elseif n == 1
        one(eltype(x)) 
    else
        ones(eltype(x), n)
    end
    
    s = State(x, fx, dfx)
    return FunIterator(N!, f, s)
end


"""
Performs a single iteration step using the update map (N!).
"""
function step!(fi::FunIterator)
    if norm(fi.S.x) > 1e12
        throw(DomainError(fi.S.x, "x is too large"))
    end
    fi.N!(fi.S)  # Update the state in-place.
end

"""
Returns the current state (x, f(x)).
"""
get_state(fi::FunIterator) = fi.S.x, fi.S.fx

"""
Sets the state variables x, fx, and dfx.  If fx or dfx are not provided, they are computed.
"""
function set_state!(fi::FunIterator, x; fx = nothing, dfx = nothing)
    fi.S.x = x

    if isnothing(fx)
        fi.S.fx = isa(fi.f, Vector) ? map(h -> h(x), fi.f) : fi.f(x)
    else
        fi.S.fx = fx
    end

    if isnothing(dfx)
        n = length(x)
        if isa(fi.f,Vector) 
            fi.S.dfx =  ones(eltype(x), n, length(fi.f))
        elseif n > 1
            fi.S.dfx = ones(eltype(x), size(fi.S.dfx))
        else
            fi.S.dfx = one(eltype(x)) 
        end
    else
        fi.S.dfx = dfx
    end

    return nothing
end

# This is where the iterations are computed until 
# the stopping criterion is met
function _get_iterations!(ds, ε, max_it)
    xn, fx = get_state(ds) 
    yy = Vector{typeof(xn)}(undef, max_it + 1)
    yy[1] = xn
    k = 1
    try
        while  norm(fx) > ε && k < max_it 
            step!(ds)
            xn, fx = get_state(ds) 
            k += 1
            yy[k] = xn
        end
    catch e
        @warn "Iteration failed at step $k: $e" 
        # @show xn
        return max_it, yy 
    end
    return k, yy
end


# Barrier function.
function stephenson_map(f::Function, g::Function, d::Int)
    if d == 1 
        return _stephenson_map(f,g)
    else 
        return _stephenson_map(f,g,d)
    end
end

# Generalized Steffenson method real values
function _stephenson_map(f::Function, g::Function)
    function N!(S::State)
        x, fx, dfx = S.x, S.fx, S.dfx
        # fx = f(x) 
        gx = g(fx) 
        fx_h = f(x + gx)
        dfx = (fx_h - fx)/gx 
        x = x - fx/dfx
        S.x = x; S.fx = f(x); S.dfx = dfx
        # return S
    end
    return N!
end


# Generalized Steffenson method for R^d → R
function _stephenson_map(f::Function, g::Function, d)
    J(x) = construct_gradient(x, f, g, d)
    function N!(S::State)
        x = S.x
        Jx, fx = J(x) 
        nJ = norm(Jx) 
        if nJ > 0 
            x_new = x - fx*Jx/nJ^2
        else 
            x_new = x
        end
        S.x = x_new; S.fx = fx
    end
    return N!
end

# Evaluate function and gradient matrix
function construct_gradient(x, f, g, d)
    J = zeros(eltype(x), d) 
    fx = f(x)
    gx = g(fx)
    G = zeros(eltype(x), d)
    for  k in 1:d 
        G .= 0.0 
        G[k] = gx
        J[k] = (f(x .+ G) - fx)/gx
    end
    return J, fx 
end

# Evaluate function and jacobian matrix
function construct_jacobian(x, fx, f, g, d)
    J = zeros(eltype(x), d, d) 
    gx = g.(fx)
    G = zeros(eltype(x), d)
 # J(x) = [ (f[n](x .+ G(x, n, k)) - f[n](x))/g(f[n](x)) for n in 1:dim, k in 1:dim] 
    for n in 1:d, k in 1:d 
        G .= 0.0 
        G[k] = gx[n]
        J[n,k] = (f[n](x .+ G) - fx[n])/gx[n]
    end
    return J 
end


# Generate a estimated jacobian function using the same technique for functions
# from R^k -> R^k._     
function stephenson_map(f::Array{Function}, g::Function, d::Int)
    J(x,fx) = construct_jacobian(x, fx, f, g, d)
    function N!(S::State)
        x = S.x; fx = S.fx
        Jx = J(x,fx) 
        if 0 < abs(det(Jx)) < Inf 
            x_new =  x - inv(Jx)*fx 
            S.x = x_new
        else
            S.x = x
        end
        S.fx = map(h -> h(S.x), f)
    end
    return N!
end

# Barrier function.
function stephenson_map_accel(f::Function, g::Function, d::Int)
    if d == 1 
        return _stephenson_map_accel(f,g)
    else 
        return _stephenson_map_accel(f,g,d)
    end
end

function _stephenson_map_accel(f::Function, g::Function)
    function N!(S::State)
        xp, fxp, dfxp = S.x, S.fx, S.dfx
        beta = -1/dfxp
        gx = g(fxp*beta)
        fx_h = f(xp + gx)
        dfx = (fx_h - fxp)/gx 
        x = xp - fxp/dfx
        # @show xp, fxp, dfx
        fx = f(x)
        S.x = x; S.fx = fx; S.dfx = dfx
        end
    return N!
end


function stephenson_map_accel_secant(f::Function, g::Function, d::Int)
    function N!(S::State)
        xp, fxp, dfxp = S.x, S.fx, S.dfx
        beta = dfxp
        gx = g(fxp*beta)
        fx_h = f(xp + gx)
        dfx = (fx_h - fxp)/gx 
        x = xp - fxp/dfx
        fx = f(x)
        beta = -(x - xp)/(fx - fxp)
        S.x = x; S.fx = fx; S.dfx = beta
        end
    return N!
end


function construct_jacobian(x, fx, Jx, f, g, d)
    J = zeros(eltype(x), d, d) 
    gx = [ g(-fx[n]/Jx[n,k]) for n in 1:d, k in 1:d ]
    G = zeros(eltype(x),d)
    # J(x) = [ (f[n](x .+ G(x, n, k)) - f[n](x))/g(f[n](x)) for n in 1:dim, k in 1:dim] 
    for n in 1:d, k in 1:d 
        G .= 0.0 
        G[k] = gx[n,k]
        J[n,k] = (f[n](x .+ G) - fx[n])/gx[n,k]
    end
    return J 
end

# Generate a estimated jacobian function using the same technique for functions
# from R^k -> R^k._     
function _stephenson_map_accel(f::Array{Function}, g::Function, d::Int)
    J(x, fx, Jx) = construct_jacobian(x, fx, Jx, f, g, d)
    function N!(S::State)
        x = S.x; fx = S.fx; Jx = S.dfx
        Jx = J(x, fx, Jx) 
        # @show x, fx, Jx
        if 0 < abs(det(Jx)) < Inf 
            x_new =  x - inv(Jx)*fx 
            S.x = x_new
        else
            S.x = x
        end
        S.fx = map(h -> h(S.x), f)
        S.dfx = Jx
    end
    return N!
end

# Generalized Steffenson method for R^d → R
function _stephenson_map_accel(f::Function, g::Function, d)
    J(x, fx, Jx) = construct_gradient(x, fx, Jx, f, g, d)
    function N!(S::State)
        x = S.x; fx = S.fx; Jx = S.dfx
        Jx = J(x, fx, Jx) 
        nJ = norm(Jx) 
        if nJ > 0 
            x_new = x - fx*Jx/nJ^2
        else 
            x_new = x
        end
        S.x = x_new; S.fx = f(x_new); S.dfx = Jx
    end
    return N!
end

# Evaluate function and gradient matrix
function construct_gradient(x, fx, Jx, f, g, d)
    J = zeros(eltype(x), d) 
    gx = g.(-fx./Jx)
    G = zeros(eltype(x), d)
    for  k in 1:d 
        G .= 0.0 
        G[k] = gx[k]
        J[k] = (f(x .+ G) - fx)/gx[k]
    end
    return J
end
