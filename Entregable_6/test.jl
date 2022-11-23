# include the module 
include("Entregable_6.jl")
using  .Entregable_6

# add function to plot the graph 
import .Entregable_6:plot_graph

# Line
V1 = 1:3; E1 = [(1, 2), (2, 3) ]
out1 = calcular_coordenadas(V1, E1)
plot_graph(out1, E1)

# Square
V2 = 1:4; E2 = [(1, 2), (2, 3), (3, 4), (4, 1)]
out2 = calcular_coordenadas(V2, E2)
plot_graph(out2, E2)

# Hexa
V3 = 1:6; E3 = [(1, 2), (1, 5), (2, 5), (2, 4), (3, 4), (3, 6), (4, 6)]
out3 = calcular_coordenadas(V3, E3)
plot_graph(out3, E3)


# load  external libraries
using Graphs, GraphPlot

g = SimpleGraph()

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


# Line
V1 = 1:3; E1 = [(1, 2), (2, 3) ]
g1 = buid_graph(V1, E1)
gplot(g1)

# Hexa
V2 = 1:6; E2 = [(1, 2), (1, 5), (2, 5), (2, 4), (3, 4), (3, 6), (4, 6)]
g2 = buid_graph(V2, E2)
gplot(g2)
