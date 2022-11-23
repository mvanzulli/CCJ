# add dependencies
using Sobol
# export 
export reachable_set

# add integrate function
include("./../Entregable_3/Ejercicio_3_1.jl")

abstract type HyperRectangle{N} end
" Rectangular box

## Fields
- `ub` upper bound point ∈ R²
- `lb` upper bound point ∈ R²

"
struct Box <: HyperRectangle{2}
    ub::Vector
    lb::Vector
    " Box constructor with center and radious"
    function Box(;c::Vector, r::Real)
        ub = c .+ r
        lb = c .- r
        new(ub,lb)
    end
    " Box constructor with center and radious"
    function Box(ub::Vector, lb::Vector)
        new(ub,lb)
    end
end


" Extracts upper point "
upper(b::Box) = b.ub 
" Extracts lower point "
lower(b::Box) = b.lb 
" Extract box dim "
dim(b::Box) = length(b.ub) 

" Computes the reachable solution set. "
function reachable_set(f::Function, alg::RK4, t₀, t₁, X0::Box; num_samples = 10_000)
    # Extract box dim 
    N = dim(X0)
    # Create SobolSeq
    s = SobolSeq(lower(X0), upper(X0))
    # Vector of vectors that contains the Y at time T for each x in s
    YT = Vector{Vector{Float64}}(undef, num_samples)
    for point in 1:num_samples 
        # compute next initial point
        x₀ᵢ = next!(s)
        # solve diff eq
        yᵢ = integrate(f, alg, t₀, t₁, x₀ᵢ)
        length(yᵢ) == 2 && throw(ArgumentError("yₜ length must be a R² vector"))        
        # save the las point
        YT[point] =yᵢ[end][1]
    end
    # build the set 
    lbₓ = minimum(getindex.(YT,1)); ubₓ = maximum(getindex.(YT,1)) 
    lbⱼ = minimum(getindex.(YT,2)); ubⱼ = maximum(getindex.(YT,2)) 

    return Box([ubₓ,ubⱼ],[lbₓ,lbⱼ])
end