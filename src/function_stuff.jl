using ForwardDiff: derivative

func_list_benchmark = [x -> (x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5),
x -> (x - 1)^3 - 1,
x ->  exp(x^2+x*cos(x)-1)*sin(x) + log(x^2+1), 
x -> abs(x^2-9)]

string_list_benchmark = [L"f_1(x) = \Pi (x - i)",
               L"f_2(x) = x -> (x - 1)^3 - 1",
               L"f_3(x) = x ->  exp(x^2+x*cos(x)-1)*sin(x) + log(x^2+1)", 
               L"f_4(x) = x -> abs(x^2-9)"]


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
    z -> 0.2*tanh(abs(z)/22)*exp(im*angle(z)), 
    z ->  min(0.2, abs(z))*exp(im*angle(z)), 
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
        x = x1[1] + im * x1[2]
        isnan(u(x)) && return SVector(real(x), imag(x))
        x_new =  x - u(x) 
        return SVector(real(x_new), imag(x_new))
    end
    return N
end

