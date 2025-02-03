using ForwardDiff: derivative
using LinearAlgebra


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

# Generalized Steffenson method 
function stephenson_map(f::Function, g::Function)
        h(z) = g(f(z))
        d(z) = (f(z + h(z)) - f(z))/h(z)
        u(z) = f(z)/d(z)
    function N(z1, p, n)
        z = z1[1] + im * z1[2]
        isnan(u(z)) && return SVector(real(z), imag(z))
        z_new =  z - u(z) 
        return SVector(real(z_new), imag(z_new))
    end
    return N
end


# Generalized Steffenson method real values
function stephenson_map_real(f::Function, g::Function)
    function N(x1)
        x = x1[1]
        fx = f(x) 
        gx = g(fx) 
        fx_h = f(x + gx)
        x_new = x - fx*gx/(fx_h - fx)
        return [x_new]
    end
    return N
end

# Evaluate function and jacobian matrix
function construct_jacobian(x,f,g,d)
    J = zeros(BigFloat, d,d) 
    fx = map(h -> h(x), f)
    G(x, n, k) = setindex!(zeros(BigFloat, d), g(fx[n]), k)  
    # J(x) = [ (f[n](x .+ G(x, n, k)) - f[n](x))/g(f[n](x)) for n in 1:dim, k in 1:dim] 
    for n in 1:d, k in 1:d 
        J[n,k] = (f[n](x .+ G(x, n, k)) - fx[n])/g(fx[n])
    end
    return J, fx 
end


# Generate a estimated jacobian function using the same technique for functions
# from R^k -> R^k._     
function stephenson_map_ndim(f::Array{Function}, g::Function, dim::Int)
    J(x) = construct_jacobian(x, f, g, dim)
    function N(x)
        Jx, fx = J(x) 
        if 0 < abs(det(Jx)) < Inf 
            x_new =  x - inv(Jx)*fx 
            return x_new
        else
            return x
        end
    end
    return N
end
