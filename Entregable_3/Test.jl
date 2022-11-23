# load  external libraries
import Test,
       DifferentialEquations as DiEq,
       Polynomials as Pol

# include the file 
include("Entregable_3.jl")
# load Entregable3 module
using .Entregable_3

###################
# Test Ejercicio 1 
###################

"Lotka Volterra diferential equation "
function lotkavolterra(x, t)
    α, β, γ, δ = 1.5, 1., 3., 1.
    [α * x[1] - β * x[1] * x[2], δ * x[1] * x[2] - γ * x[2]]
end

# problem def
x₀ = ones(2); t₀ = .0; t₁ = .5

Xt = integrate(lotkavolterra, RK4(h=0.001), t₀, t₁, x₀)

"Lotka Volterra for Differential Equations solution "
function lotkavolterra_diffeq!(du, u, p, t)
    α, β, γ, δ = 1.5, 1., 3., 1.
    du[1] = α * u[1] - β * u[1] * u[2]
    du[2] = δ * u[1] * u[2] - γ * u[2]
end

# create problem
prob = DiEq.ODEProblem(lotkavolterra_diffeq!, x₀, (t₀, t₁))
# solve it 
sol_difeq = DiEq.solve(prob, DiEq.RK4())

Test.@testset "Ejercicio 1" begin
        Test.@test sol_difeq.u[end] ≈ Xt[end-1][1] atol = 1e-2
end

###################
# Test Ejercicio 2 
###################

Test.@testset "Ejercicio 2" begin 
  
    # Test b basis  
    pout1 = bernstein_basis(3, 1) 
    p_test1 = Pol.Polynomial([0, 3, -6, 3])
    diff = sum(abs(pout1(x) - p_test1(x)) for x in rand(1_000))
    Test.@test diff ≈ 0 atol = 1e-10
   
    # Test b coefficients
    p_test2 = Pol.Polynomial([3, 2, -5])
    Test.@test bernstein_coefficients(p_test2) == [3, 4, 0]
end


