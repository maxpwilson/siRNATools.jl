
struct RNA_Alphabet
    bases::Array{Char, 1}
    bits::Array{UInt8, 1}
    base_to_bit::Dict{Char, UInt8}
    bit_to_base::Dict{UInt8, Char}
    RNA_Alphabet(bases) = new(bases, [UInt8(x) for x in 0:length(bases) - 1], Dict(zip(bases, [UInt8(x) for x in 0:length(bases) - 1])), Dict(zip([UInt8(x) for x in 0:length(bases) - 1], bases)))
end
Base.length(T::RNA_Alphabet) = length(T.bases)
function Base.iterate(T::RNA_Alphabet , (el, i) = (T.bases[1], 1))
    if i > length(T) 
        return nothing
    elseif i == length(T)
        return (T.bases[i], (nothing, i + 1))
    else
        return (T.bases[i], (T.bases[i+1] , i+1))
    end
end


const RNA_ALPHABET = RNA_Alphabet(['A', 'C', 'G', 'U', 'N'])
const BASES = RNA_ALPHABET.base_to_bit
const BIT_BASES = RNA_ALPHABET.bit_to_base
