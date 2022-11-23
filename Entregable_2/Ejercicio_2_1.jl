# Structs defs
struct ADN
    dat::String
    function ADN(dat::String)
        # check data is only formed by nuc
        neucleotidos = ['A', 'C', 'G', 'T']
        char_founded = 0
        for n in neucleotidos
            char_founded += length(findall(n, dat))
        end
        @assert char_founded == length(dat) "data must be formed only by $neucleotidos"
        new(dat)
    end
end
struct ARN
    dat::String
    function ARN(dat::String)
        # check data is only formed by nuc
        neucleotidos = ['A', 'C', 'G', 'U']
        char_founded = 0
        for n in neucleotidos
            char_founded += length(findall(n, dat))
        end
        @assert char_founded == length(dat) "data must be formed only by $neucleotidos"
        new(dat)
    end        
end
# transcribir fucntion
function transcribir(adn::ADN)
    uracil = 'U'
    timina = 'T'
    ARN(reverse(replace(adn.dat, timina => uracil)))
end

