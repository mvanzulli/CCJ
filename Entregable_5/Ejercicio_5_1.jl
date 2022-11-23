import LinearAlgebra: dot, norm

export Vector1E, val, val_index

"Vector with a single not null entry:

## Fields:
- `vᵢ` value 
- `i` index value 
- `N` vector length "
struct Vector1E
    vᵢ::Number
    i::Int
    n::Int
    Vector1E(x; i, n) = 0 < i ≤ n && new(x, i, n)
    # TODO: Add show method into the struct
end
function Base.show(io::IO, v::Vector1E)
    print(io, "1 index vector: = \n 0 ∀ i ≠ $(val_index(v)) \n $(val(v)) if i == $(val_index(v)) ")
end
"Returns the value of the vector `v`"
val(v::Vector1E) = v.vᵢ
"Returns index of the non-zero value of the vector `v`"
val_index(v::Vector1E) = v.i  

Base.length(v::Vector1E) = v.n
Base.getindex(v::Vector1E, idx::Int) = idx == val_index(v) ? val(v) : 0 
Base.getindex(v::Vector1E, idxs::Vector{Int}) =  [ getindex(v,i) for i in idxs ] 
Base.iterate(v::Vector1E, state=1) =  (1 > state || state > length(v)) ? nothing : (v[state], state + 1)
Base.size(v::Vector1E) = (v.n,)
Base.size(v::Vector1E, dim) = dim == 1 ? length(v) : 1
Base.:-(v::Vector1E) = Vector1E(-val(v); i = val_index(v), n = length(v))

function Base.:+(v1::Vector1E, v2::Vector1E)
    @assert length(v1) == length(v2) throw(DimensionMismatch("dim v1 must be the same as v2"))
    if val_index(v1) == val_index(v2) 
        return Vector1E(val(v1) + val(v2), i = val_index(v1), n = length(v1))    
    else
        v = zeros(length(v))
        v[val_index(v1)] = val(v1) 
        v[val_index(v2)] = val(v2) 
        return v
    end
end

function Base.:+(v1::Vector{T}, v2::Vector1E) where T<:Real
    @assert length(v1) == length(v2) throw(DimensionMismatch("dim v1 must be the same as v2"))
    v1 = copy(V1)
    v1[val_index(v2)] += val(v2)
    return v1
end

function Base.:+(v1::Vector1E, v2::Vector{T}) where T<:Real
    @assert length(v1) == length(v2) throw(DimensionMismatch("dim v1 must be the same as v2"))
    v2 = copy(v2)
    v2[val_index(v1)] += val(v1)
    return v2
end

function Base.:-(v1::Vector{T}, v2::Vector1E) where T<:Real
    @assert length(v1) == length(v2) throw(DimensionMismatch("dim v1 must be the same as v2"))
    v1 = copy(v1)
    v1[val_index(v2)] -= val(v2)
    return v1
end

function Base.:-(v1::Vector1E, v2::Vector{T}) where T<:Real
    @assert length(v1) == length(v2) throw(DimensionMismatch("dim v1 must be the same as v2"))
    v2 = copy(v2)
    v2[val_index(v1)] -= val(v1)
    return -v2
end

function Base.:*(m::Matrix{T}, v::Vector1E) where T <:Number
    @assert size(m, 2) == length(v) throw(DimensionMismatch("cols of M must be row of v"))
    return val(v) * m[:, val_index(v)]
end

# Overlead LinearAlgebra methods
dot(v1::Vector1E, v2::Vector{T}) where T<:Real = v2[val_index(v1)] * val(v1)
norm(v1::Vector1E) where T<:Real = abs(val(v1))

