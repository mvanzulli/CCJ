# Obj/Func exported
export integrate, RK4
#=
Implementación de Range-Kutta o4 (x' = f(x(t),t))

xₖ₊₁ = h ₛ∑ₖ₌₁ bₖwₖ

donde:
bₖ ∈ (1/6,2/6,2/6,1/6)
y
w ∈
w₁ = f(x,t) 
w₂ = f(x + (h * w₁ / 2, t + h / 2))
w₃ = f(x + (h * w₂ / 2, t + h / 2))
w₄ = f(x + h * w₃, t + h) 

=#
" Range Kutta oreder 4 struct

## Fields
-`h` time step
-`b` weights
"
Base.@kwdef struct RK4
    h::Real
    b::NTuple{4, Rational} = (1/6, 2/6, 2/6, 1/6) 
end

" Devuelve el paso `h` de Rnage-Kutte-4"
step(rk4::RK4) = rk4.h

" Devuelve el vector `b` de Rnage-Kutte-4"
coefs(rk4::RK4) = rk4.b

" Devuelve el vector w "
function compute_wᵥ(xₖ, f::Function, rk4::RK4, t::Real)
    # extract step 
    h = step(rk4)
    # weights
    w₁ = f(xₖ,t) 
    w₂ = f(xₖ .+ h / 2 .* w₁ , t + h / 2) 
    w₃ = f(xₖ .+ h / 2 .* w₂, t + h / 2)
    w₄ = f(xₖ .+ h .* w₃, t + h) 
    wᵥ = (w₁, w₂, w₃, w₄)    
end

"Devuelve el valor de xₖ₊₁"
function compute_xₖ₊₁(xₖ, f::Function, rk4::RK4, t)
    
    # Extraer parámetros del método 
    h = step(rk4)
    bᵥ = coefs(rk4)

    # compute w vector
    wᵥ = compute_wᵥ(xₖ, f, rk4, t)
    
    # inicializar xₖ₊₁
    xₖ₊₁ = xₖ
    # compute xₖ₊₁
    for (i,bₖ) in enumerate(bᵥ)
        wₖ = wᵥ[i] 
        xₖ₊₁ += h * bₖ * wₖ
    end
    return xₖ₊₁
end

"Integra el método numérico entre t₀ y t₁"
function integrate(f::Function, rk4::RK4, t₀::Real, t₁::T, x₀) where T<:Real
    # inicialización 
    t = t₀; iₜ = 1 
    xₖ = x₀; hist_x = Vector{Vector{typeof(x₀)}}([])
    global t, hist_x, xₖ
    # integración
    while t ≤ t₁
        # calculo el siguiente paso k+1
        xₖ₊₁ = compute_xₖ₊₁(xₖ, f, rk4, t)
        # gruardo resultado
        push!(hist_x, [xₖ₊₁]) 
        # increment and update
        t += step(rk4)
        xₖ = xₖ₊₁
    end
    return hist_x
end
