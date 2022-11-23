module Entregable2 
    # cedula de identidad
    const CI = 46718847
    # Ejercicio_2_1
    include("Ejercicio_2_1.jl")
    export transcribir
    
    # Ejercicio_2_2
    include("Ejercicio_2_2.jl")
    export traducir

    # Ejercicio_2_3
    include("Ejercicio_2_3.jl")
    export estimar_pi

end

using .Entregable2, Test
@testset "Entregable 2" begin
    # Test 2_1
    adn = Main.Entregable2.ADN("CCTAGGACCAGGTT")
    arn = Main.Entregable2.ARN("UUGGACCAGGAUCC")
    @test transcribir(adn) == arn 
    
    # Test_2_2
    out = traducir(Main.Entregable2.ARN("CCU"))
    @test out.dat == [Main.Entregable2.Aminoacido("Pro")]
    # Test_3_3
    buf = Main.Entregable2.Buffon(10, 1)
    πestim = estimar_pi(buf,100009)
    @test 2.9 ≤ πestim ≤ 3.6
end


