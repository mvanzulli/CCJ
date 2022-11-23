# Obj/Func exported
export bernstein_basis,
       bernstein_coefficients,
       bernstein_enclosure

# include dependencies
import DynamicPolynomials as DP
import Polynomials as Pol

" Computes the Berenstein polynomial bases Bᵢ "
function bernstein_basis(l::Int, i::Int)
    # declare variable x
    DP.@polyvar x
    i<0 || i> l ? Bᵢ = 0 : Bᵢ = binomial(l,i) * x^i * (1-x)^(l-i)
    # convert 
    Bᵢ = convert(DP.Polynomial, Bᵢ)    
end

" Computes Berenstein coefficients bᵢ in the base Bᵢ  "
function bernstein_coefficients(p::Pol.Polynomial)
    # convert 
    p = convert(DP.Polynomial, p)    
    # p coefs
    a = reverse(DP.coefficients(p))
    # l and i 
    l = DP.maxdegree(p); i = l;
    #
    b = Vector{Rational}(undef, l + 1)
    for i in 0:l
        bᵢ = .0
        for j in 0:i
            bᵢ+= binomial(i,j) / binomial(l,j) * a[j+1]
        end
        b[i+1] = bᵢ
    end
    return b
end

# add domain generalization
function bernstein_coefficients(p::Pol.Polynomial, X::NTuple{2,Real})
    # convert 
    p = convert(DP.Polynomial, p)    
    # p coefs
    a = reverse(DP.coefficients(p))
    # l and i 
    l = maxdegree(p); i = l;
    #
    b = Vector{Rational}(undef, l + 1)
    for i in 0:l
        bᵢ = .0
        for j in 0:i
            for k in j:l
                bᵢ+= binomial(i,j) / binomial(l,j) * (X[2] - X[1])^j *
                     binomial(j,k) * X[1]^(k-j) * a[k+1]
            end
        end
        b[i+1] = bᵢ
    end
    return b
end
" Computes polynomial enclosure given by the Bernstein Polynomials for a polynomial `p` in the `X`` domain"
function bernstein_enclosure(p::Pol.Polynomial, X::Tuple{Number,Number})
    # convert 
    p = convert(DP.Polynomial, p)
    # check X domain 
    @assert X[1] ≤ X[2] "X₁ must be larger than X₂ but are X₁ = $X₁ and X₂ = $X₂"
    # compute b coefs
    b = bernstein_coefficients_new(p, X) 
    # return max and min 
    return extrema(b)
end 

" Function to convert `pol` form Polynomial to DynamicPolynomials "
function Base.:convert(::Type{DP.Polynomial}, pol::Pol.Polynomial)
    DP.@polyvar x
    sum(x^(i-1) * p for (i, p) in enumerate(pol.coeffs))
end