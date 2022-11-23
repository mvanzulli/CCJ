# Problem sketch
                    #H
# y = D    # ------------------------------    #- 
                       # / L                   #- D
# y = 0    # -----------------------------     #-

# Buffon struct
Base.@kwdef struct Buffon
    D::Float64=10
    L::Float64=8
end 

linespace(buf::Buffon) = buf.D
nead_length(buf::Buffon) = buf.L

# Point 
struct Point
    x::Real
    y::Real
end
# Horizontal Line
struct HorizontalLine
    y::Real
end
yline = yval(hl::HorizontalLine) = hl.y 
# Overload + 
function Base.:+(p1::Point, p2::Point)
    Point(p1.x + p2.x, p1.y + p2.y)
end 
# Overload * for a scalar 
function Base.:*(α::Number, p::Point)
    Point(α*p.x, α*p.y)
end 

# Neadle 
struct Neadle 
    start::Point
    finish::Point
end
start(n::Neadle) = n.start
finish(n::Neadle) = n.finish
# Overload maximum for each neadle coord 
function Base.:maximum(n::Neadle)
    # xₘₐₓ
    xₘₐₓ = maximum((start(n).x, finish(n).x))
    # yₘₐₓ
    yₘₐₓ = maximum((start(n).y, finish(n).y))
    return (xₘₐₓ, yₘₐₓ)
end 
# Overload minimum for each neadle coord 
function Base.:minimum(n::Neadle)
    # xₘᵢₙ
    xₘᵢₙ = minimum((start(n).x, finish(n).x))
    # yₘᵢₙ
    yₘᵢₙ = minimum((start(n).y, finish(n).y))
    return (xₘᵢₙ, yₘᵢₙ)
end 
# overload cap function
function Base.:∩(n::Neadle, hl::HorizontalLine)
    # extract y line 
    y_line = yval(hl)
    # Compute max and min y of the neadle
    yₘₐₓ = maximum(n)[2]
    yₘᵢₙ = minimum(n)[2]
    return yₘᵢₙ ≤ y_line ≤ yₘₐₓ 
end
# overload length function
function Base.:length(n::Neadle) # in order to not imoprt LinearAlgebra norm function length is employed
    sn = start(n)
    fn = finish(n)
    sqrt((fn.x - sn.x)^2 + (fn.y - sn.y)^2)
end

"""Funcion to estimate π"""
function estimar_pi(buf::Buffon, N=1000)
    # Create horizontal lines
    yD = HorizontalLine(linespace(buf))
    y0 = HorizontalLine(0)
    # Generate slope samples 
    p = 2π*rand(N)
    # Horizontal distance where the samples can live, however this has no effect since y = const  
    H = 2.0        
    # Generate origin points 
    starts = [Point(xᵢ, yᵢ) for (xᵢ, yᵢ) in zip(H*rand(N), linespace(buf)*rand(N)) ]
    # Generate relative finshes points 
    local_finish =  [nead_length(buf)*Point(cos(pᵢ), sin(pᵢ)) for pᵢ in p]
    # Build finshes
    finishes =  starts .+ local_finish
    # Create neadles vector
    samples = [Neadle(startᵢ, finishᵢ) for (startᵢ, finishᵢ) in zip(starts, finishes)]
    # Check if the neadles cross the hlines
    neadles_crossed = 0
    for neadle in samples 
        # check length
        !(length(neadle) ≈ nead_length(buf)) && "The nealde length is not as it is expected" 
        # return true if yD is interescted
        crossyD = ∩(neadle, yD)
        # return true if y0 is interescted
        crossy0 = ∩(neadle, y0)
        # check only one is crossed
        crossyD && crossy0 && throw(ArgumentError("L must  ≤ D"))
        # add to the counter 
        crossyD || crossy0 ? neadles_crossed +=1 : continue
    end
    # sample freq 
    f = neadles_crossed / N
    # 
    πestim = 2*nead_length(buf)/ (f*linespace(buf))
    return πestim

end