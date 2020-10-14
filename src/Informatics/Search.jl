"""
    calculate_Peq(::String) :: Array{UInt64, 1}

Takes as input an RNA strand and translates it into an array of Unsigned Integers each one representing one letter and each bit referring to a position.
Used as a pre-processing step for Myers searching algorithm.  Maximum length of searched pattern is 64.

### Example

AACGCU becomes [110000, 001010, 000100, 000001]
"""
function calculate_Peq(pattern::String) :: Array{UInt64, 1}
    m = length(pattern)
    Peq :: Array{UInt64, 1} = zeros(UInt64, length(RNA_ALPHABET))
    for i in 1:m
        y = pattern[i]
        for x in RNA_ALPHABET
            if x == y
                Peq[BASES[x] + 1] |= one(UInt64) << (i-1)
            end
        end
    end
    Peq
end

"""
    motif_to_transcript_match(Peq::Array{Unit64, 1}, m::Int64, refseq::ReferenceSequence, min_K::Int) :: Array{Tuple{UInt64, UInt64}}

Takes as input Peq calculated in [calculate\\_Peq](@ref), m = the length of the pattern to search for, the ReferenceSequence to search in, and min\\_k =
the maximum distance to output.  Implements the string matching algorithm described in _A Fast Bit-Vector Algorithm for Approximate String Matching Based
on Dynamic Programming_ by Gene Myers.  Output is a tuple of all Levenshtein distances less than min_k and their positions in refseq.
"""
function motif_to_transcript_match(Peq::Array{UInt64, 1}, m::Int64,  refseq::ReferenceSequence, min_k::Int) :: Array{Tuple{Int, Int}}
    min_j = 0
    out::Array{Tuple{Int, Int}} = []
    n::UInt64 = refseq.length
    Pv::UInt64 = (one(UInt64) << m) - one(UInt64)
    Mv::UInt64 = zero(UInt64)
    Xv::UInt64 = 0
    Mh::UInt64 = 0
    Ph::UInt64 = 0
    Xh::UInt64 = 0
    Eq::UInt64 = 0
    dist = m
    for j in 1:n
        @inbounds Eq = Peq[get_refseq_pos(refseq, j) + 1]
        Xv = Eq | Mv
        Xh = (((Eq & Pv) + Pv) ⊻ Pv) | Eq

        Ph = Mv | ~(Xh | Pv)
        Mh = Pv & Xh
        if (Ph >> (m-1)) & 1 != 0
            dist += 1
        elseif (Mh >> (m-1)) & 1 != 0
            dist -= 1
        end
        if dist < min_k
            push!(out, (dist, j))
        end

        Ph <<= 1
        Mh <<= 1
        Pv = Mh | ~(Xv | Ph)
        Mv = Ph & Xv
    end
    out
end

"""
    find_match_sequences(motif::String, sequence::String, mismatches::Int) :: Array{String, 1}

Function takes as input a motif, sequence, and number of mismatches to search for.  Output is an array of all substrings of sequence which have a
Hamming distance of exactly mismatches to motif.
"""
function find_match_sequences(motif::String, sequence::String, mismatches::Int, tail::Int=0, head::Int=0) :: Array{String, 1}
    mtchs::Array{String, 1} = []
    for i in 1:(length(sequence) - length(motif) + 1)
        (i - tail > 0) ? start = i - tail : start = 1
        (i + length(motif) - 1 + head <= length(sequence)) ? stop = i + length(motif) - 1 + head : stop = length(sequence)
        (HammingDist(motif, sequence[i:i+length(motif) - 1]) == mismatches) && push!(mtchs, sequence[start:stop])
    end
    return mtchs
end


"""
    mismatch_positions(::String, ::String) :: Array{Int, 1}

Function returns an array of all positions in which strings used as input differ.  Strings must be the same length.
"""
function mismatch_positions(seq1::String, seq2::String) :: Array{Int, 1}
    out::Array{Int, 1} = []
    temptotal::Array{Int, 1} = []
    if length(seq1) == length(seq2)
        for i in 1:length(seq1)
            (seq1[i] != seq2[i]) && push!(out, i)
        end
    elseif length(seq1) > length(seq2)
        temptotal = [x for x in 1:length(seq2)]
        out = [x for x in 1:length(seq1) - length(seq2) + 1]
        for j in 1:(length(seq1) - length(seq2) + 1)
            temp::Array{Int, 1} = []
            for i in j:length(seq2) + j - 1
                (seq1[i] != seq2[i - j + 1]) && push!(temp, i)
            end
            (length(temp) < length(temptotal)) && (temptotal = copy(temp))
        end
        out = vcat(out, temptotal)
    else
        temptotal = [x for x in 1:length(seq1)]
        out = [x for x in 1:length(seq2) - length(seq1) + 1]
        for j in 1:(length(seq2) - length(seq1) + 1)
            temp::Array{Int, 1} = []
            for i in j:length(seq1) + j - 1
                (seq1[i - j + 1] != seq2[i]) && push!(temp, i)
            end
            (length(temp) < length(temptotal)) && (temptotal = copy(temp))
        end
        out = vcat(out, temptotal)
    end
    return out
end
