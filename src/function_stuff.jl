using ForwardDiff: derivative


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
        h(x) = g(f(x))
        d(x) = (f(x + h(x)) - f(x))/h(x)
        u(x) = f(x)/d(x)
    function N(x1, p, n)
        x = x1[1]
        isinf(abs(x)) && return SVector(x)
        isnan(u(x)) && return SVector(x)
        x_new =  x - u(x) 
        return SVector(x_new)
    end
    return N
end


# Generate a estimated jacobian function using the same technique for functions
# from R^k -> R^k._     
function stephenson_map_ndim(f::Array{Function}, g::Function, dim::Int)
    G(x, n, k) = setindex!(zeros(dim), g(f[n](x)), k)  
    J(x) = [ (f[n](x .+ G(x, n, k)) - f[n](x))/g(f[n](x)) for n in 1:dim, k in 1:dim] 
    function N(x, p, n)
        if any(isinf.(abs.(x)))
            return SVector{dim}(x)
        end
        fx = map(h -> h(x), f)

        if any(isinf.(abs.(fx)))
            return SVector{dim}(x)
        end
        Jx = J(x) 
        if any(isnan.(Jx)) 
            Jx[isnan.(Jx)] .= 1.
        end

        # @show Jx, fx, x
        if 0 < abs(det(Jx)) < Inf 
            x_new =  x - inv(Jx)*fx 
            return SVector{dim}(x_new)
        else
            return SVector{dim}(x)
        end
    end
    return N
end
