
export CustomSet

"A custom set struct with a vector of elements:

## Fields:
- `elements` value "
struct CustomSet
    elements::Vector{T} where T<:Real
end

" Returns the elements of an specific `cs` CustomSet "
elements(cs::CustomSet) = cs.elements
Base.length(cs::CustomSet) = length(elements(cs))
Base.getindex(cs::CustomSet,i::Int) = 1 ≤ i ≤ length(cs) ? elements(cs)[i] : throw(BoundsError("i is out of range"))

Base.:(==)(cs1::CustomSet, cs2::CustomSet ) = unique(elements(cs1)) == unique(elements(cs2))
for f in (:(Base.:∩),:(Base.:∪)) 
    @eval $f(cs1::CustomSet, cs2::CustomSet) = CustomSet($f(elements(cs1), elements(cs2)))
end
for f in (:(Base.:∈),:(Base.:∉))
    @eval $f(x::Number, cs::CustomSet) = $f(x, elements(cs))
end
Base.push!(cs::CustomSet, x::Number) = x ∉ elements(cs) ? CustomSet(push!(elements(cs), x)) : cs
Base.sum(cs::CustomSet) = sum(elements(cs))
Base.iterate(cs::CustomSet, step = 1) = 1 ≤ step ≤ length(cs) ? (elements(cs)[step], step + 1) : nothing