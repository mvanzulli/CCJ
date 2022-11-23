# load  external libraries
using Test, LinearAlgebra

# include the module 
include("Entregable_5.jl")
using  .Entregable_5

Test.@testset "Ejercicio_5_1" begin 

    v1E = Vector1E(6, i=2, n=10)
    @test v1E == Vector1E(6, i=2, n=10)
    @test length(v1E) == 10 
    @test -v1E == Vector1E(-6, i=2, n=10)
    @test v1E + v1E == Vector1E(12, i=2, n=10)
    @test ones(10, 10 ) * v1E  == val(v1E)*ones(10)
    @test norm(v1E) == 6
    @test dot(v1E, v1E) == norm(v1E)^2
end

Test.@testset "Ejercicio_5_2" begin 
    @test CustomSet([1, 2, 3, 1]) == CustomSet([1, 2, 3])
    @test CustomSet([1, 2, 3]) ∩ CustomSet([1, 2]) == CustomSet([1, 2])
    @test CustomSet([1, 2, 3]) ∪ CustomSet([5, 6]) == CustomSet([1, 2, 3, 5, 6])
    @test 1 ∈ CustomSet([1, 2, 3])
    @test -1 ∉ CustomSet([1, 2, 3])
    @test length(CustomSet([1, 2, 3])) == 3
    @test push!(CustomSet([1, 2, 3]), 5) == CustomSet([1, 2, 3, 5])
    @test push!(CustomSet([1, 2, 3]), 1) == CustomSet([1, 2, 3])
    @test hasmethod(iterate, Tuple{CustomSet})
    @test sum(CustomSet([1, 2, 3])) == 6
end
