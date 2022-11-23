    
# single character to genereate chains
chars = ['U', 'C', 'A', 'G']
# chars permutations
perms = [c1 * c2 * c3 for c1 in chars for c2 in chars for c3 in chars ]

# indexes in perms who starts with 
indexes_first_U = startswith.(perms, chars[1])
indexes_first_C = startswith.(perms, chars[2])
indexes_first_A = startswith.(perms, chars[3])
indexes_first_G = startswith.(perms, chars[4])

# inline function to check the second character of a spring vector
secondis(s::String, char::Char) = s[2] == char

# indexes in perms who's second letter is   
indexes_second_U = secondis.(perms, chars[1])
indexes_second_C = secondis.(perms, chars[2])
indexes_second_A = secondis.(perms, chars[3])
indexes_second_G = secondis.(perms, chars[4])

# indexes in perms who's final letter is 
indexes_ends_U = endswith.(perms, chars[1])
indexes_ends_C = endswith.(perms, chars[2])
indexes_ends_A = endswith.(perms, chars[3])
indexes_ends_G = endswith.(perms, chars[4])

# Indexes for each aminoacid in perms
indexes_Phe = indexes_first_U .* indexes_second_U .* (indexes_ends_U .+ indexes_ends_C ) 
indexes_Leu = indexes_first_C .* indexes_second_U .+ indexes_first_U .* indexes_second_U .* (indexes_ends_A .+ indexes_ends_G ) 
indexes_Ser = indexes_first_U .* indexes_second_C .+
                indexes_first_A .* indexes_second_G .*(indexes_ends_U .+ indexes_ends_C) 
indexes_Tyr = indexes_first_U .* indexes_second_A
indexes_Cys = indexes_first_U .* indexes_second_G .* (indexes_ends_U .+ indexes_ends_C)
indexes_Trp = indexes_first_U .* indexes_second_G .* indexes_ends_G 
indexes_Pro = indexes_first_C .* indexes_second_C 
indexes_His = indexes_first_C .* indexes_second_A .* (indexes_ends_U .+ indexes_ends_C)
indexes_Gin = indexes_first_C .* indexes_second_A .* (indexes_ends_A .+ indexes_ends_G)
indexes_Arg = indexes_first_C .* indexes_second_G .+ indexes_first_A .* indexes_second_G .*(indexes_ends_A .+ indexes_ends_G)
indexes_Lle = indexes_first_A .* indexes_second_U .* (indexes_ends_U .+ indexes_ends_C .+ indexes_ends_A) 
indexes_Met = indexes_first_A .* indexes_second_U .* indexes_ends_G 
indexes_Thr = indexes_first_A .* indexes_second_C 
indexes_Asn = indexes_first_A .* indexes_second_G .* (indexes_ends_U .+ indexes_ends_C)
indexes_Lis = indexes_first_A .* indexes_second_G .* (indexes_ends_A .+ indexes_ends_G)
indexes_Val = indexes_first_G .* indexes_second_U 
indexes_Ala = indexes_first_G .* indexes_second_C 
indexes_Asp = indexes_first_G .* indexes_second_A .* (indexes_ends_U .+ indexes_ends_C)
indexes_Glu = indexes_first_G .* indexes_second_A .* (indexes_ends_A .+ indexes_ends_G)
indexes_stop = indexes_first_U .* indexes_second_A .* (indexes_ends_A .+ indexes_ends_G) +
                indexes_first_U .* indexes_second_G .* (indexes_ends_A)

# Chains
amin_table = Dict{String, Vector{String}}( )
chains_Phe = perms[findall(indexes_Phe.==1)]; amin_table["Phe"] = chains_Phe
chains_Leu = perms[findall(indexes_Leu.==1)]; amin_table["Leu"] = chains_Leu 
chains_Ser = perms[findall(indexes_Ser.==1)]; amin_table["Ser"] = chains_Ser
chains_Tyr = perms[findall(indexes_Tyr.==1)]; amin_table["Tyr"] = chains_Tyr
chains_Cys = perms[findall(indexes_Cys.==1)]; amin_table["Cys"] = chains_Cys
chains_Trp = perms[findall(indexes_Trp.==1)]; amin_table["Trp"] = chains_Trp
chains_Pro = perms[findall(indexes_Pro.==1)]; amin_table["Pro"] = chains_Pro
chains_His = perms[findall(indexes_His.==1)]; amin_table["His"] = chains_His
chains_Gin = perms[findall(indexes_Gin.==1)]; amin_table["Gin"] = chains_Gin
chains_Arg = perms[findall(indexes_Arg.==1)]; amin_table["Arg"] = chains_Arg
chains_Lle = perms[findall(indexes_Lle.==1)]; amin_table["Lle"] = chains_Lle
chains_Met = perms[findall(indexes_Met.==1)]; amin_table["Met"] = chains_Met
chains_Thr = perms[findall(indexes_Thr.==1)]; amin_table["Thr"] = chains_Thr
chains_Asn = perms[findall(indexes_Asn.==1)]; amin_table["Asn"] = chains_Asn
chains_Lis = perms[findall(indexes_Lis.==1)]; amin_table["Lis"] = chains_Lis
chains_Val = perms[findall(indexes_Val.==1)]; amin_table["Val"] = chains_Val
chains_Ala = perms[findall(indexes_Ala.==1)]; amin_table["Ala"] = chains_Ala
chains_Asp = perms[findall(indexes_Asp.==1)]; amin_table["Asp"] = chains_Asp
chains_Glu = perms[findall(indexes_Glu.==1)]; amin_table["Glu"] = chains_Glu
chains_stop = perms[findall(indexes_stop.==1)]; amin_table["stop"] = chains_stop

# Amino struct
struct Aminoacido
    chain::String
    function Aminoacido(chain) 
        length(chain) > 3 && throw(ArgumentError("An aminoacid must contain less than 3 letters and contains $(length(chain))")) 
        new(chain)
    end
end

# Protain chain struct
struct CadenaProteica
    dat::Vector{Aminoacido}
end

function traducir(arn::ARN, tab = amin_table )
    # extract data
    dat = arn.dat
    # output vector 
    out = CadenaProteica([])
    # number of amin in the ARN
    num_amin = div(length(arn.dat), 3)
    # check length is 3
    length(dat) % 3 !== 0 && throw(ArgumentError("The arn must contain a length multiple of 3 but is $(length(dat))")) 
    # split dat in 3 char strings
    chains = [dat[3(i -1) + 1:3(i)] for i in 1:num_amin] 
    for chain in chains
        for (amin, vec_chains) in tab
            if chain in vec_chains
                push!(out.dat, Aminoacido(amin))
                break
            end
        end
    end                
    return out         
end