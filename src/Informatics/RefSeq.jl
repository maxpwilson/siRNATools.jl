
"""
Structure for reference sequences.  Compresses RNA data into 2 bits of information from the 8 of a normal character string.  Can only use bases A, C, G, U.
- A => 00
- C => 01
- G => 10
- U => 11
"""
struct ReferenceSequence
    data::Vector{UInt64}
    nmask::BitArray{1}
    length::UInt64
end


"""
    encode_refseq(::String) :: ReferenceSequence

Function that takes a string containing only A, C, G, and U and encodes it into a ReferenceSequence type
"""
function encode_refseq(seq)
    data = Vector{UInt64}(undef, cld(length(seq), 32))
    nmask = falses(length(seq))
    i = 1
    for j in 1:lastindex(data)
        x = UInt64(0)
        r = 0
        while r < 64 && i <= lastindex(seq)
            nt = seq[i]
            if nt in ['A', 'C', 'G', 'U']
                x |= convert(UInt64, BASES[nt]) << r
            else
                nmask[i] = true
            end
            i += 1
            r += 2
        end
        data[j] = x
    end
    ReferenceSequence(data, nmask, length(seq))
end


"""
    decode_refseq(::ReferenceSequence) :: String

Function takes a ReferenceSequence type and returns a String containing the RNA bases.
"""
function decode_refseq(refseq::ReferenceSequence)
    out = ""
    i = 1
    for d in refseq.data
        for r in 0:2:62
            if refseq.nmask[i] 
                out *= 'N'
            else
                out *= BIT_BASES[d >> r & 0b11]
            end
            (i==refseq.length) && return out
            i += 1
        end
    end
    out
end

"""
    decode_refseq_partial(::ReferenceSequence, ::UnitRange) :: String

Function which takes a ReferenceSequence type and outputs specificied range of bases as a String.
"""
function decode_refseq_partial(refseq::ReferenceSequence, rg::UnitRange)::String
    out = ""
    for j in rg
        out *= BIT_BASES[get_refseq_pos(refseq, j)]
    end
    out
end


"""
    get_refseq_pos(::ReferenceSequence, ::Int)

Function to retrieve bit value of base at specificied position in a reference sequence.
"""
function get_refseq_pos(refseq, pos)
    @assert(pos >= 1 && pos <= refseq.length)
    i = cld(pos, 32)
    j = (pos - (32 * (i - 1)) - 1) * 2
    refseq.data[i] >> j & 0b11
end
