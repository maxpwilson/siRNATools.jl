module Specificity
using CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter, BSON, Base.Threads
using BSON: @save, @load

include("Path.jl")
include("RNAAlphabet.jl")
include("RefSeq.jl")
include("Maintenance.jl")


"""
    reverse_complement(::String) :: String

Takes the reverse complement of an RNA string
"""
function reverse_complement(pattern::String) :: String
    pattern = uppercase(pattern)
    BASES = ['A', 'C', 'G', 'U']
    @assert sum([x in BASES for x in pattern]) == length(pattern)
    complements = ['U', 'G', 'C', 'A']
    b_to_c = Dict{Char, Char}(zip(BASES, complements))
    r_c = ""
    for base in pattern
        r_c = b_to_c[base] * r_c
    end
    return r_c
end

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
function motif_to_transcript_match(Peq::Array{UInt64, 1}, m::Int64,  refseq::ReferenceSequence, min_k::Int) :: Array{Tuple{UInt64, UInt64}}
    min_j = 0
    out::Array{Tuple{UInt64, UInt64}} = []
    n = refseq.length
    Pv::UInt64 = (one(UInt64) << m) - one(UInt64)
    Mv::UInt64 = zero(UInt64)
    dist = m
    for j in 1:n
        Eq = Peq[get_refseq_pos(refseq, j) + 1]
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
function find_match_sequences(motif::String, sequence::String, mismatches::Int) :: Array{String, 1}
    mtchs::Array{String, 1} = []
    for i in 1:(length(sequence) - length(motif) + 1)
        (evaluate(Hamming(), motif, sequence[i:i+length(motif) - 1]) == mismatches) && push!(mtchs, sequence[i:i+length(motif) - 1])
    end
    return mtchs
end

"""
    mismatch_positions(::String, ::String) :: Array{Int, 1}

Function returns an array of all positions in which strings used as input differ.  Strings must be the same length.
"""
function mismatch_positions(seq1::String, seq2::String) :: Array{Int, 1}
    out::Array{Int, 1} = []
    for i in 1:length(seq1)
        if length(seq2) < i 
            push!(out, i)
        else
            (seq1[i] != seq2[i]) && push!(out, i)
        end
    end
    return out
end


"""
    find_genome_matches(pattern::String, excluded_gene::String = "", verbose::Bool = true, minimum_matches = 5) :: Array{Tuple{String, Int64}}

Function takes as input a pattern to search the genome for, excluded_gene to exclude the gene of interest from the search, verbose set to true displays
a progress bar showing progress of search, and minimum\\_matches which is the amount of mismatches searched for + 1.  Output is an array of tuples of transcript names
and the lowest Hamming distance found to the pattern within that transcript.
"""
function find_genome_matches(pattern::String, excluded_gene::String = "",  verbose::Bool = true, minimum_matches = 5, pbar_string::String = "Searching Genome... ") :: Array{Tuple{String, Int64}}
    (length(ALLREFSEQ) == 0) && return []
    out::Array{Tuple{String, Int64}} = []
    (verbose ==true) && (p = Progress(length(ALLREFSEQ), 0.1, pbar_string))
    for (name, T) in ALLREFSEQ
        (excluded_gene != "") && ((name in GENETRANSCRIPTS[excluded_gene]) && continue)
        min_match = minimum_matches
        min_pos = 0
        matches = motif_to_transcript_match(calculate_Peq(pattern),length(pattern), T, minimum_matches)
        for (dist, pos) in matches
            start = (pos-length(pattern)+1 > 0) ? (pos-length(pattern)+1) : 1
            stop = (pos <= T.length) ? pos : T.length
            match = decode_refseq_partial(T, start:stop)
            if length(mismatch_positions(pattern, match)) < min_match
                min_match = length(mismatch_positions(pattern, match))
                min_pos = pos
            end
        end
        (min_match < minimum_matches) && push!(out, (name, min_match))
        (verbose == true) && ProgressMeter.next!(p)
    end
    out
end

"""
    compress_genome_matches(::Array{Tuple{String, Int64}}) :: Dict{String, Array{Int64, 1}}

Function takes as input the output of [find\\_genome\\_matches](@ref) and collapses it into a dictionary of Gene names and an array of all lowest Hamming mismatch 
of the transcripts of that gene.
"""
function compress_genome_matches(raw_data::Array{Tuple{String, Int64}}) :: Dict{String, Array{Int64, 1}}
    out = Dict{String, Array{Int64, 1}}()
    (length(TRANSCRIPTGENE) == 0) && return out
    for (name, match) in raw_data
        if !(haskey(out, TRANSCRIPTGENE[name]))
            out[TRANSCRIPTGENE[name]] = [match]
        else
            push!(out[TRANSCRIPTGENE[name]], match)
        end
    end
    return out
end


"""
    final_calc(pattern::String, raw_data::Array{Tuple{String, Int64}}, compressed_data::Dict{String, Array{Int64, 1}})

Function takes as input the pattern being searched for, the raw\\_data from [find\\_genome\\_matches](@ref), and the compressed data from [compress\\_genome\\_matches](@ref).
Ouput is a Tuple containing the mismatch_counts dictionary which adds up the number of genes with a minimum Hamming distance of 0-4, and the specificity score.
"""
function final_calc(pattern::String, raw_data::Array{Tuple{String, Int64}}, compressed_data::Dict{String,Array{Int64, 1}})
    mismatch_counts = Dict{Int64, Int64}([x => 0 for x in 0:4])
    gene_correction = Dict()
    transcript_list = []
    min_score::Float64 = 5
    min_match::Int64 = 5
    for (gene, matches) in compressed_data
        filter!(x -> x == minimum(matches), matches)
        mismatch_counts[minimum(matches)] += 1
        (minimum(matches) < min_match) && (min_match = minimum(matches))
        gene_correction[gene] = []
    end
    for (name, match) in raw_data
        if minimum(match) == min_match
            match_patterns = find_match_sequences(pattern, decode_refseq(ALLREFSEQ[name]), min_match)
            for match_pattern in match_patterns
                score::Float64 = 0
                for mismatch in mismatch_positions(pattern, match_pattern)
                    if mismatch in 1:7
                        score += 1
                    elseif mismatch in 8:9
                        score += 1.2
                    elseif mismatch == 10
                        score += 1.25
                    elseif mismatch == 11
                        score += 1.5
                    else
                        score += 1.9
                    end
                end
                (score < min_score) && (min_score = score)
            end
        end
    end
    #for (gene, corrections) in gene_correction
    #    (length(corrections) > 0) && (mismatch_counts[length(corrections[1])] += (length(unique(corrections)) - 1))
    #end
    (min_score == 2.9) && (min_score = 3)
    (mismatch_counts, min_score)
end

"""
    Calculate_Specificity(patterns, excluded_gene="", rg=2:18, verbose=true) :: DataFrame

Function takes as input the pattern to be searched against the genome, the excluded gene to be ignored, the range of pattern being used, and a boolean which
controls whether progress bars will be shown.  Output is a DataFrame with a column for the pattern, the number of genes with minimum mismatch distance of 0-4
and the specificity score.
"""
function Calculate_Specificity(patterns, excluded_gene::String="", rg::UnitRange{Int64} = 2:18, verbose::Bool=true) :: DataFrame
    string_array = Array{String, 1}()
    for pattern in patterns
        push!(string_array, pattern)
    end
    Calculate_Specificity(string_array, excluded_gene, rg, verbose)
end
function Calculate_Specificity(patterns::Array{String, 1}, excluded_gene::String="", rg::UnitRange{Int64} = 2:18, verbose::Bool=true) :: DataFrame
    df = DataFrame(ID=Int[], Pattern=String[], Zero=Int64[], One=Int64[], Two=Int64[], Three=Int64[], Four=Int64[], Score=Float64[])
    counter = 0
    (verbose == true) && ((p = Progress(length(patterns), 0.2, "Calculating Specificity ... ")))
    (verbose == true) && (update!(p, 0))
    @threads for i in 1:length(patterns)
        pattern = patterns[i]
        counter += 1
        RP = reverse_complement(pattern[rg])
        raw_data = find_genome_matches(RP, excluded_gene, false, 5, "Searching strand $(counter) of $(length(patterns)) ... ")
        compressed_data = compress_genome_matches(raw_data)
        (mismatchs, spec_score) = final_calc(RP, raw_data, compressed_data)
        push!(df, [i, pattern, mismatchs[0], mismatchs[1], mismatchs[2], mismatchs[3], mismatchs[4], spec_score])
        (verbose == true) && ProgressMeter.next!(p)
    end
    sort!(df, :ID)
    df
end

export Calculate_Specificity, ReferenceSequence

end
