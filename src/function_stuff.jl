using ForwardDiff: derivative
using LinearAlgebra

"""
    mutable struct FunIterator{N <: Number}

A mutable struct that facilitates iteration of a function. It stores the function, the current state (`x`),
and the function value at the current state (`fx`). The type parameter `N` specifies the numeric type of the state.

# Fields
- `N::Function`: The function to be iterated.  It should take `x` and return a tuple of `(new_x, f(new_x))`.
- `x::Union{Vector{N}, N}`: The current state of the system (scalar or vector).
- `fx::Union{Vector{N}, N}`: The function value at the current state (scalar or vector).
"""
mutable struct State{N <: Number}
    x::Union{Vector{N}, N}
    fx::Union{Vector{N}, N}
    dfx::Union{Vector{N}, N}
end

mutable struct FunIterator
    N!::Function
    S::State
end


function FunIterator(f::Function, x) 
    s = State(x, x, x)
    return FunIterator(f, s)
end
    

function step!(fi::FunIterator)
    if norm(fi.S.x) > 1e5
        throw(DomainError(fi.S.x, "x is too big"))
    end
     fi.N!(fi.S)
end

function get_state(fi::FunIterator)
    return fi.S.x, fi.S.fx
end

function set_state!(fi::FunIterator, x) 
        fi.S.x = x
end

function set_state!(fi::FunIterator, x, fx) 
        fi.S.x = x
        fi.S.fx = fx
end
# This is where the iterations are computed until 
# the stopping criterion is met
function _get_iterations!(ds, ε, max_it)
    xn_1, fx = get_state(ds) 
    yy = Vector{typeof(xn_1)}(undef, max_it + 1)
    yy[1] = xn_1
    step!(ds)
    xn, fx = get_state(ds) 
    k = 1
    # stopping criterion is ∥x_n - x_{n-1}∥ + ∥f(x_{n-1})∥ ≤ ε
    while norm(xn - xn_1) + norm(fx) > ε  
        (k > max_it) && break 
        yy[k+1] = xn
        xn_1 = xn
        try 
            step!(ds)
        catch 
            k = max_it + 1 
            break 
        end
        xn, fx = get_state(ds) 
        k += 1
    end
    return k, yy
end

function ∂f(f)
# Warning. This trick works only for holomorphic functions. 
    ∂f∂z(z) = (g(x) = f(x+z); derivative(x->real(g(x)),0) 
            + im * derivative(x->imag(g(x)),0))
    return ∂f∂z
end

function N_map(z, f, ∂f∂z)
    dz = f(z)/∂f∂z(z)
    return  z - dz
end

# For modified two step newton method with β
function beta_map(f, β)
    ∂f∂z = ∂f(f)
    N(z) = N_map(z, f, ∂f∂z)
    function N_β(z1, p, n)
        z = z1[1] + im * z1[2]
        N_z = N(z)
        z_new =  N_z - β * f(N_z)/∂f∂z(z)
        return SVector(real(z_new), imag(z_new))
    end
    return N_β
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
        fx = f(x) 
        gx = g(fx) 
        fx_h = f(x + gx)
        dfx = (fx_h - fx)/gx 
        x = x - fx/dfx
        S.x = x; S.fx = fx; S.dfx = dfx
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
    G(k) = setindex!(zeros(eltype(x), d), gx, k)  
    for  k in 1:d 
        J[k] = (f(x .+ G(k)) - fx)/gx
    end
    return J, fx 
end

# Evaluate function and jacobian matrix
function construct_jacobian(x,f,g,d)
    J = zeros(eltype(x), d, d) 
    fx = map(h -> h(x), f)
    gx = g.(fx)
    G(n, k) = setindex!(zeros(eltype(x), d), gx[n], k)  
    # J(x) = [ (f[n](x .+ G(x, n, k)) - f[n](x))/g(f[n](x)) for n in 1:dim, k in 1:dim] 
    for n in 1:d, k in 1:d 
        J[n,k] = (f[n](x .+ G(n, k)) - fx[n])/gx[n]
    end
    return J, fx 
end


# Generate a estimated jacobian function using the same technique for functions
# from R^k -> R^k._     
function stephenson_map(f::Array{Function}, g::Function, d::Int)
    J(x) = construct_jacobian(x, f, g, d)
    function N!(S::State)
        x = S.x
        Jx, fx = J(x) 
        if 0 < abs(det(Jx)) < Inf 
            x_new =  x - inv(Jx)*fx 
            S.x = x_new
        else
            S.x = x
        end
        S.fx = fx
    end
    return N!
end

