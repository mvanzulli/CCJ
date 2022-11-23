# Declare vertex_coords type  
const Tuple2Real = NTuple{2,T} where {T<:Real}

# Use external dependencies  
using LazySets: Singleton, LineSegment
using Graphs: SimpleGraph, floyd_warshall_shortest_paths, add_edge!
using LinearAlgebra: norm, dot
using Plots

# Export requested functions
export calcular_coordenadas


" Function build graph with a set of E and V"
function buid_graph(V, E)
    # Initialize
    nᵥ = length(V)
    g = SimpleGraph(nᵥ)

    # add edges
    for edge in E
        add_edge!(g, edge...)
    end

    # return the builded graph
    return g
end

" Return dᵢⱼ  "
compute_dᵢⱼ(g::SimpleGraph) = floyd_warshall_shortest_paths(g).dists

" Return lᵢⱼ  "
compute_lᵢⱼ(dᵢⱼ, L) = dᵢⱼ * L

" Return dᵢⱼ  "
compute_kᵢⱼ(dᵢⱼ, K) = (dᵢⱼ .^ -2) * K

" Computes the number of vertex"
nvertex(V) = length(V)

"Computes the vector `Δ`"
compute_Δ(vertex_coords, kᵢⱼ, lᵢⱼ) = [norm(compute_∇E(i, vertex_coords, kᵢⱼ, lᵢⱼ)) for i in 1:length(vertex_coords)]

" Updates the vector `Δ` "
function update_Δ!(Δ, vertex_coords, kᵢⱼ, lᵢⱼ)
    for i in 1:length(vertex_coords)
        Δ[i] = norm(compute_∇E(i, vertex_coords, kᵢⱼ, lᵢⱼ))
    end
end


" Initialize the algorithm "
function initialize(V, E, L=1, K=1)

    g = buid_graph(V, E)
    dᵢⱼ = compute_dᵢⱼ(g)
    kᵢⱼ = compute_kᵢⱼ(dᵢⱼ, L)
    lᵢⱼ = compute_lᵢⱼ(dᵢⱼ, K)

    # initialize vertex_coords
    vertex_coords = Vector{Tuple2Real}(undef, nvertex(V))
    [vertex_coords[vertex] = Tuple(rand(2)) for vertex in 1:nvertex(V)]
    # [vertex_coords[vertex] = (vertex, vertex + 1) for vertex in 1:nvertex(V)]

    return vertex_coords, dᵢⱼ, kᵢⱼ, lᵢⱼ, 0
end

"Computes the energy gradient `∇Eₘ`"
function compute_∇E(m::Integer, vertex_coords::Vector{Tuple2Real}, kᵢⱼ::Matrix{R}, lᵢⱼ::Matrix{Q}) where {R<:Real,Q<:Real}
    # Initialize ∇E vector
    ∇Eₘ = zeros(Real, 2)
    # Extract coordinates of vertex m
    Xₘ = vertex_coords[m]
    # Compute derivatives ∂E∂X
    for i in 1:length(vertex_coords)

        if i != m # add zero to the gradient since the sum is for i ≆ m 

            # compute difᵢ = (xₘ - xᵢ) and  difⱼ = (yₘ - yᵢ)
            Xᵢ = vertex_coords[i]
            
            
            dif = Xₘ .- Xᵢ
            ∇Eₘ .+= kᵢⱼ[m, i] .* (dif .- lᵢⱼ[m, i] .* dif ./ norm(dif))

        end

    end
    return ∇Eₘ
end

" Compte the energy Hessian `HE` "
function compute_HE(m::Integer, vertex_coords::Vector{Tuple2Real}, kᵢⱼ::Matrix{R}, lᵢⱼ::Matrix{Q}) where {R<:Real,Q<:Real}
    
    # Initialize ∇E vector
    HEₘ = zeros(Real, (2, 2))
    
    # Extract coordinates of vertex m
    Xₘ = vertex_coords[m]
    
    # Compute derivatives ∂E∂X
    for i in 1:length(vertex_coords)

        if i == m

            HEₘ .+= zeros(2, 2)

        else

            # compute difᵢ = (xₘ - xᵢ) and  difⱼ = (yₘ - yᵢ)
            Xᵢ = vertex_coords[i]
            dif = Xₘ .- Xᵢ
            HEₘ[1, 1] += kᵢⱼ[m, i] * (1 - lᵢⱼ[m, i] * (dif[2])^2 / (dot(dif, dif))^(3 / 2))
            HEₘ[2, 2] += kᵢⱼ[m, i] * (1 - lᵢⱼ[m, i] * (dif[1])^2 / (dot(dif, dif))^(3 / 2))
            diagonal = kᵢⱼ[m, i] * lᵢⱼ[m, i] * (dif[1] * dif[2]) / (dot(dif, dif))^(3 / 2)
            HEₘ[1, 2] += diagonal
            HEₘ[2, 1] += diagonal

        end

    end

    return HEₘ

end

" Computes the graph new coordinates. "
function calcular_coordenadas(
    V::UnitRange,
    E::Vector{NTuple{2,T}},
    ϵ=1e-8,
    tol_iters=400) where {T<:Integer}

    # Initialize graph 
    L = 1; K = 1
    vertex_coords, dᵢⱼ, kᵢⱼ, lᵢⱼ, iter = initialize(V, E)

    # Find maximum index of Δ at the first step
    Δ = compute_Δ(vertex_coords, kᵢⱼ, lᵢⱼ)
    Δₘ, idx_Δmax = findmax(Δ)

    # Execute the algorithm
    while Δₘ ≥ ϵ
        # Set the max to the particle p        
        Δₚ = Δₘ

        while Δₚ ≥ ϵ && iter ≤ tol_iters
        
            # Initialize ∇E and HE
            ∇E = compute_∇E(idx_Δmax, vertex_coords, kᵢⱼ, lᵢⱼ)
            HE = compute_HE(idx_Δmax, vertex_coords, kᵢⱼ, lᵢⱼ)
        
            # compute increment 
            δ = HE \ -∇E
            
            # update postion 
            vertex_coords[idx_Δmax] =  Tuple(vertex_coords[idx_Δmax] .+ δ)
        
            # update Δ, Δₚ, and iter 
            update_Δ!(Δ, vertex_coords, kᵢⱼ, lᵢⱼ)
            Δₚ = Δ[idx_Δmax]
            iter += 1
            
        end
        
        # Select the new max
        Δₘ, idx_Δmax = findmax(Δ)

    end

    return vertex_coords

end

" Plot a graph given coordinates in a vector of tuples `vertex_coords` "
function plot_graph(vertex_coords, conectivity)
    # select backend
    fig = plot()
    gr()
    theme(:juno)
    plot_font = "Computer Modern"
    default(
        fontfamily=plot_font,
        linewidth=2,
        framestyle=:box,
        label=nothing,
        grid=true,
    )
    # plot nodes
    for vertex in vertex_coords
        fig = plot!(Singleton(vertex...), seriestype=:scatter, markersize= 12)
    end
    # plot line segments
    for conec in conectivity
        fig = plot!(
            LineSegment(
                vertex_coords[conec[1]] |> collect,
                vertex_coords[conec[2]]|> collect,
                    ),
                )
    end

    # deploy
    display(fig)
end
