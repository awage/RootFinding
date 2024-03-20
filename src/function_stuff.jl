using ForwardDiff: derivative

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
x -> exp((x-3)^2 + 7(x - 3) - 30) - 1,
x -> (x + sin(2/x) * x*x), 
x -> (x-1/4)^(1/3) ]

string_list = [L"f(x) = (x^2 - 1)(x^2 + 1)",
L"f(x) = x^3 - 1",
L"f(x) = x^{12} - 1",
L"f(x) = (x^2 - 4)(x + 1.5)(x - 0.5)",
L"f(x) =(x+2)(x+1.5)^2 (x-0.5)(x-2)",
L"f(x) = \sin(x)",
L"f(x) = (x - 1)^3 + 4(x-1)^2 - 10",
L"f(x) = \sin(x-14/10)^2 - (x - 14/10)^2 + 1",
L"f(x) = x^2 - e^x - 3x + 2",
L"f(x) = \cos(x-3/4) - x + 3/4",
L"f(x) = (x + 1)^3 - 1",
L"f(x) = (x-2)^3 - 10",
L"f(x) = (x + 5/4)~  e^{(x + 5/4)^2} - \sin(x + 5/4)^2 + 3\cos(x + 5/4) + 5",
L"f(x) = exp((x-3)^2 + 7(x - 3) - 30) - 1",
L"f(x) = (x + sin(2/x)  x^2)", 
L"f(x) = (x - 1/4)^{1/3} "]


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

function beta_map(f)
    ∂f∂z = ∂f(f)
    N(z) = N_map(z, f, ∂f∂z)
    function N_β(z1, p, n)
        β = p[1]
        z = z1[1] + im * z1[2]
        N_z = N(z)
        z_new =  N_z - β * f(N_z)/∂f∂z(z)
        return SVector(real(z_new), imag(z_new))
    end
    return N_β
end


function beta_map_df(f)
    ∂f∂z = ∂f(f)
    N(z) = N_map(z, f, ∂f∂z)
    function N_β(z1, p, n)
        β = p[1]
        z = z1[1] + im * z1[2]
        N_z = N(z)
        z_new =  N_z - β(∂f∂z(z), N_z - z) * f(N_z)/∂f∂z(z)
        return SVector(real(z_new), imag(z_new))
    end
    return N_β
end

function beta_map_real(f)
    ∂f∂x(x) = derivative(f,x)
    N(x) = N_map(x, f, ∂f∂x)
    function N_β!(dx, x, p, n)
        β = p[1]
        Nx = N(x[1])
        dx[1] = Nx - β*f(Nx)/∂f∂x(x[1])
        return
    end
    return N_β!
end
