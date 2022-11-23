# Add libraries
using StatsPlots, Distributions
# Exprots
export simular

struct Accion # Aₜ
    distribucion::Normal
end

struct Juego # [A₁ A₂ A₃ ]
    acciones::Vector{Accion}
end


function Juego(; k::Int = 10, μ = zeros(k), σ = ones(k))
    vec = [Accion(Normal(μi, σi)) for (μi, σi) in zip(μ, σ)]
    Juego(vec)
end

function jugar(m::Juego, i::Int)
    p = m.acciones[i]
    rand(p.distribucion)
end

Base.@kwdef struct UCB
    c::Float64 = 2
end

function simular(J::Juego, alg::UCB; budget::Int = 1000)
    # Numero de acciones
    k = length(J.acciones)

    # Inicializar contadores.
    Q = zeros(k) # vector de recompensas 
    N = zeros(k) # número de veces que la accion A fue seleccionada
    acciones = Vector{Int}(undef,0)
    recompensas_promedio = Vector{Float64}(undef,0)
    t = 1 
    # parameter c UCB
    c = alg.c

    while t ≤ budget
        # First explore 
        
        A = if minimum(N) == 0
            idx = findall(iszero, N)
            rand(idx)       
        else
            argmax(Q + c * ( N / log(t)).^(-.5))
        end
        push!(acciones, A)
        R = jugar(J, A)

        # Actualizar contadores.
        N[A] += 1
        Q[A] += 1 / N[A] * (R - Q[A])
        push!(recompensas_promedio, Q[A])
        t += 1
    end
    return acciones, recompensas_promedio
end

