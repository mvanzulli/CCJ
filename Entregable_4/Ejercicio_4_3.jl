
# Import dependencies
import Polynomials as Pol
import DynamicPolynomials as DP
using BenchmarkTools

# Obj/Func exported
export bernstein_basis, bernstein_coefficients, bernstein_enclosure

" Function to convert `pol` form Polynomial to DynamicPolynomials "
function Base.:convert(::Type{DP.Polynomial}, pol::Pol.Polynomial)
    DP.@polyvar x
    sum(x^(i-1) * p for (i, p) in enumerate(pol.coeffs))
end

#############
# Naive way
#############
"Naive method for Bernstein coefficients computation"
struct Naive end

" Computes Berenstein coefficients bᵢ in the base Bᵢ  "
function bernstein_coefficients(p::Pol.Polynomial, ::Naive)::Vector
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
function bernstein_coefficients(p::Pol.Polynomial, X::NTuple{2,Real}, ::Naive)::Vector
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
function bernstein_enclosure(p::Pol.Polynomial, X::Tuple{Number,Number}, alg::Naive)
    # convert 
    p = convert(DP.Polynomial, p)
    # check X domain 
    @assert X[1] ≤ X[2] "X₁ must be larger than X₂ but are X₁ = $X₁ and X₂ = $X₂"
    # compute b coefs
    b = bernstein_coefficients(p, X, Naive()) 
    # return max and min 
    return extrema(b)
end 

"Fast method for Bernstein coefficients computation"
struct Fast end


function bernstein_coefficients(p::Pol.Polynomial, ::Fast)::Vector
    # convert 
    p = convert(DP.Polynomial, p)    
    # p coefs
    a = reverse(DP.coefficients(p))
    # p max degree
    l = DP.maxdegree(p);
    # max degree of each monomial 
    maxmono = DP.maxdegree.(reverse(DP.monomials(p)))
    
    # brenstein coefs vector
    b = zeros(Rational, l + 1)
    
    # Fill the last l + 1
    b[l + 1] = a[l + 1]
    pop!(maxmono)
  
    # compute coefs and sum
   for k in maxmono
        for i in k : l
            b[i + 1] += binomial(l - k, i - k) / binomial(l, i) * a[k + 1] 
        end
    end
    return b
end


function bernstein_coefficients(p::Pol.Polynomial, X::NTuple{2,Real}, ::Fast)::Vector
    # convert 
    p = convert(DP.Polynomial, p)     
    a = reverse(DP.coefficients(convert(DP.Polynomial, p)))
    # p max degree
    l = DP.maxdegree(p);
    # max degree of each monomial 
    maxmono = DP.maxdegree.(reverse(DP.monomials(p)))

    # brenstein coefs vector
    b = zeros(Rational, l + 1)

   # compute and sum
   for k in maxmono
        for i in 0:l
            if k == l
                b[i+1] += X[1]^(k-i) * X[2]^i * a[i + 1]
            else
                for j in max(0,k - l - k):min(k, i)
                    b[i+1] += binomial(l - k, i - j) * binomial(k, j) / binomial(l, i) *
                    X[1]^(k - j) * X[2]^j * a[j + 1] # a[j +1] why it doesn't work? 
                end
            end
        end
    end
    return b
end

function bernstein_enclosure(p::Pol.Polynomial, X::Tuple{Number,Number}, ::Fast)::NTuple{2,Number}
    # convert 
    p = convert(DP.Polynomial, p)
    # check X domain 
    @assert X[1] ≤ X[2] "X₁ must be larger than X₂ but are X₁ = $X₁ and X₂ = $X₂"
    # compute b coefs
    b = bernstein_coefficients(p, X, ::Fast()) 
    # return max and min 
    return extrema(b)
end


