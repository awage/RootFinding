using ForwardDiff: derivative

func_list_benchmark = [x -> (x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5),
x -> (x - 1)^3 - 1,
x ->  exp(x^2+x*cos(x)-1)*sin(x) + log(x^2+1), 
x -> abs(x^2-9),
x -> abs(x^2-9)]

string_list_benchmark = [L"f_1(x) = \Pi (x - i)",
               L"f_2(x) = x -> (x - 1)^3 - 1",
               L"f_3(x) = x ->  exp(x^2+x*cos(x)-1)*sin(x) + log(x^2+1)", 
               L"f_4(x) = x -> abs(x^2-9)",
               L"f_4(x) = x -> abs(x^2-9)"]

root_list_benchmark = [2. , 2. , 0. , 3, -3 ] 
x0_list_benchmark = [1.5, 1.5, 0.35, 2.8, -10]

func_list = [x -> (x*x - 1) * (x*x + 1),
x -> x*x*x - 1,
x -> x^12 - 1,
x -> (x*x - 4)*(x + 1.5)*(x - 0.5),
x -> (x+2)*(x+1.5)^2*(x-0.5)*(x-2),
x -> sin(x),
x -> (x - 1)^3 + 4 * (x-1)^2 - 10,
x -> sin(x-14/10)^2 - (x - 14/10)^2 + 1,
x -> x*x - exp(x) - 3x + 2,
x -> cos(x-3/4) - x + 3/4,
x -> (x + 1)^3 - 1,
x -> (x-2)^3 - 10,
x -> (x + 5/4) * exp((x + 5/4)*(x + 5/4)) - sin((x + 5/4))^2 + 3 * cos((x + 5/4)) + 5,
x -> (x + sin(2/x) * x*x), 
]

string_list = [L"f_1(x) = (x^2 - 1)(x^2 + 1)",
L"f_2(x) = x^3 - 1",
L"f_3(x) = x^{12} - 1",
L"f_4(x) = (x^2 - 4)(x + 1.5)(x - 0.5)",
L"f_5(x) =(x+2)(x+1.5)^2 (x-0.5)(x-2)",
L"f_6(x) = \sin(x)",
L"f_7(x) = (x - 1)^3 + 4(x-1)^2 - 10",
L"f_8(x) = \sin(x-14/10)^2 - (x - 14/10)^2 + 1",
L"f_9(x) = x^2 - e^x - 3x + 2",
L"f_{10}(x) = \cos(x-3/4) - x + 3/4",
L"f_{11}(x) = (x + 1)^3 - 1",
L"f_{12}(x) = (x-2)^3 - 10",
L"f_{13}(x) = (x + 5/4)~  e^{(x + 5/4)^2} - \sin(x + 5/4)^2 + 3\cos(x + 5/4) + 5",
L"f_{14}(x) = (x + sin(2/x)  x^2)"]

    fam_list = [
    z -> 0.2*tanh(abs(z)/0.2)*exp(im*angle(z)), 
    z ->  min(0.2, abs(z))*exp(im*angle(z)), 
    z -> z] 

    fam_list_real = [
    z -> 0.1*tanh(z/0.1), 
    z -> min(0.1, z), 
    z -> z] 

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
        # β = p[1]
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
        # x = x1[1] + im * x1[2]
        x = x1[1]
        isinf(x) && return SVector(x)
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
        if any(isinf.(x))
            return SVector{dim}(x)
        end
        fx = map(h -> h(x), f)
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



F_list =[ x -> x^3 - 9x^2 + 28x - 30, 
    x -> sin(x) + x * cos(x), 
    x -> exp(x^2) - exp(sqrt(2) * x), 
    x -> (sin(x) - 0.5)^2, 
    x -> cos(x) - x, 
    x -> atan(x), 
    x -> (x - 1)^5 - 1, 
    x -> 4 * sin(x) - x + 1,
    x -> (x*x - 1) * (x*x + 1),
    x -> x*x*x - 1,
    x -> x^12 - 1,
    x -> (x*x - 4)*(x + 1.5)*(x - 0.5),
    x -> (x+2)*(x+1.5)^2*(x-0.5)*(x-2),
    x -> sin(x),
    x -> (x - 1)^3 + 4 * (x-1)^2 - 10,
    x -> sin(x-14/10)^2 - (x - 14/10)^2 + 1,
    x -> x*x - exp(x) - 3x + 2,
    x -> cos(x-3/4) - x + 3/4,
    x -> (x + 1)^3 - 1,
    x -> (x-2)^3 - 10,
    x -> (x + 5/4) * exp((x + 5/4)*(x + 5/4)) - sin((x + 5/4))^2 + 3 * cos((x + 5/4)) + 5,
    x -> x + sin(2/x) * x*x, ] 

F_X0 = [1, 0.5, 1, 3, 0.6, 0.5, 1.7, 2.5, 
       2. , 2., 2., 2.,
       2. , 2., 2., 2.,
       2. , 2., 2., 2.,
       2. , 2.
      ] 

F2_list = [[ x->  x[1] + exp(x[2]) - cos(x[2]), x->  3x[1] - x[2] - sin(x[2]) ], 
[x->   exp(x[1]^2) + 8x[1] * sin(x[2]), x->  x[1] + x[2] - 1 ],
[x-> sin(x[1]) + x[2] * cos(x[1]), x-> x[1] - x[2] ],
[ x-> x[1]^2 - 2x[1] - x[2] + 0.5, x-> x[1]^2 + 4x[2]^2 - 4.0 ],
[ x-> exp(x[1]^2) - exp(sqrt(2) * x[1]), x-> x[1] - x[2] ],
[x -> x[2] * x[3] + x[4] * (x[2] + x[3]),
x -> x[1] * x[3] + x[4] * (x[1] + x[3]),
x -> x[1] * x[2] + x[4] * (x[1] + x[2]),
x -> x[1] * x[2] + x[1] * x[3] + x[2] * x[3] - 1]]


F2_X0 = [[1,1], [0.2, 0.8], [0.4, 0.4], [0.5, 0.5], [-1, -3], [0.6, 1.6, 0.6, -0.2]]

 # F(x) = x^3 - 9x^2 + 28x - 30, \quad \alpha = 3 \text{ simple root.} \\
 # F(x) = \sin(x) + x \cos(x), \quad \alpha = 0 \text{ simple zero.} \\
 #  F(x) = \exp(x^2) - \exp\left(\sqrt{2}x\right), \quad \alpha = 0 \text{ simple zero.} \\
 #  F(x) = \left(\sin(x) - \frac{1}{2}\right)^2, \quad \alpha = 0 \text{ double zero.} \\
 #  F(x) = \cos(x) - x, \quad \alpha = 0.73908513 \text{ simple zero.} \\
 #  F(x) = \arctan(x), \quad \alpha = 0 \text{ simple zero.} \\
 #    \text{(g)} \quad & F(x) = (x - 1)^5 - 1, \quad \alpha = 2 \text{ simple root.} \\
 #  F(x) = 4 \sin(x) - x + 1, \quad \text{whose zeros are } \alpha_1 = -2.21008394, \alpha_2 = -0.34218505 \text{ and } \alpha_3 = 2.70206137. \\
 #  F(x_1, x_2) = \begin{pmatrix} x_1 + \exp(x_2) - \cos(x_2) \\ 3x_1 - x_2 - \sin(x_2) \end{pmatrix}, \quad \alpha = (0, 0)^T. \\
 #  F(x_1, x_2) = \begin{pmatrix} \exp(x_1^2) + 8x_1 \sin(x_2) \\ x_1 + x_2 - 1 \end{pmatrix}, \quad \alpha_1 = (0.175599, 0.824401)^T, \alpha_2 = (0.704247, 0.295753)^T. \\
 #  F(x_1, x_2) = \begin{pmatrix} \sin(x_1) + x_2 \cos(x_1) \\ x_1 - x_2 \end{pmatrix}, \quad \alpha = (0, 0)^T. \\
 #    \text{(l)} \quad & F(x_1, x_2) = \begin{pmatrix} x_1^2 - 2x_1 - x_2 + 0.5 \\ x_1^2 + 4x_2^2 - 4 \end{pmatrix}, \quad \alpha = (-0.222215, 0.993808)^T. \\
 #    \text{(m)} \quad & F(x_1, x_2) = \begin{pmatrix} \exp(x_1^2) - \exp(\sqrt{2}x_1) \\ x_1 - x_2 \end{pmatrix}, \quad \alpha = (0, 0)^T. \\
 #  F(x_1, x_2) = \begin{pmatrix} x_1 + \exp(x_2) - \cos(x_2) \\ 3x_1 - x_2 - \sin(x_2) \end{pmatrix}, \quad \alpha = (0, 0)^T. \\
 #  F(x) = (f_1(x), f_2(x), \ldots, f_n(x)), \quad \text{where } x = (x_1, x_2, \ldots, x_n)^T \text{ and } f_i : \mathbb{R}^n \to \mathbb{R}, i = 1, 2, \ldots, n, \text{ such that} \\
 #    & f_i(x) = x_i x_{i+1} - 1, \quad i = 1, 2, \ldots, n-1, \\
 #    & f_n(x) = x_n x_1 - 1. \\
 #    & \text{When } n \text{ is odd, the exact zeros of } F(x) \text{ are } \alpha_1 = (1, 1, \ldots, 1) \text{ and } \alpha_2 = (-1, -1, \ldots, -1). \text{ Results appearing in Table 2 are obtained for } n = 1001. \\
 #    \text{(p)} \quad & F(x) = (f_1(x), f_2(x), f_3(x), f_4(x)), \quad \text{where } x = (x_1, x_2, x_3, x_4)^T \text{ and } f_i : \mathbb{R}^4 \to \mathbb{R}, i = 1, 2, 3, 4, \text{ such that} \\
 #    & f_1(x) = x_2 x_3 + x_4 (x_2 + x_3), \\
 #    & f_2(x) = x_1 x_3 + x_4 (x_1 + x_3), \\
 #    & f_3(x) = x_1 x_2 + x_4 (x_1 + x_2), \\
 #    & f_4(x) = x_1 x_2 + x_1 x_3 + x_2 x_3 - 1.
# \end{align}

