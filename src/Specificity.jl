using CSV, DataFrames, StatsBase, StringDistances
using BSON: @save, @load

PATH = "C:\\Users\\mwilson\\Notebooks\\Specificity\\"
ALLT = Dict{String, String}()
GENETRANSCRIPTS = Dict{String, Array{String, 1}}()
TRANSCRIPTGENE = Dict{String, String}()


function load_RefSeq(path::String=PATH)
    if length(ALLT) == 0 
        @load "$path\\Human_mRNA_allT.bson" allT
        for (k, v) in allT
            ALLT[k] = v
        end
    end
    @assert length(allT) == length(ALLT)
    if length(GENETRANSCRIPTS) == 0
        @load "$path\\Human_mRNA_GeneTranscripts.bson" GeneTranscripts
        for (k, v) in GeneTranscripts
            GENETRANSCRIPTS[k] = v
        end
    end
    @assert length(GENETRANSCRIPTS) == length(GeneTranscripts)
    if length(TRANSCRIPTGENE) == 0
        @load "$path\\Human_mRNA_TranscriptGene.bson" TranscriptGene
        for (k, v) in TranscriptGene
            TRANSCRIPTGENE[k] = v
        end
    end
    @assert length(TranscriptGene) == length(TRANSCRIPTGENE)
end

function unload_RefSeq()
    for (k, v) in ALLT
        delete!(ALLT, k)
    end
    @assert length(ALLT) == 0
    for (k, v) in GENETRANSCRIPTS
        delete!(GENETRANSCRIPTS, k)
    end
    @assert length(GENETRANSCRIPTS) == 0
    for (k, v) in TRANSCRIPTGENE
        delete!(TRANSCRIPTGENE, k)
    end
    @assert length(TRANSCRIPTGENE) == 0
end

function reverse_complement(pattern::String) :: String
    pattern = uppercase(pattern)
    bases = ['A', 'C', 'G', 'U']
    @assert sum([x in bases for x in pattern]) == length(pattern)
    complements = ['U', 'G', 'C', 'A']
    b_to_c = Dict{Char, Char}(zip(bases, complements))
    r_c = ""
    for base in pattern
        r_c = b_to_c[base] * r_c
    end
    return r_c
end

function motif_to_transcript_match(motif::String, sequence::String) :: UInt8
    mtch::UInt8 = length(motif) + 1
    for i in 1:(length(sequence) - length(motif) + 1)
        new_mtch::UInt8 = evaluate(Hamming(), motif, sequence[i:i+length(motif) - 1])
        (new_mtch < mtch) && (mtch = new_mtch)
        (mtch == 0) && return 0
    end
    return mtch
end

function find_match_sequences(motif::String, sequence::String, mismatches::Int) :: Array{String, 1}
    mtchs::Array{String, 1} = []
    for i in 1:(length(sequence) - length(motif) + 1)
        (evaluate(Hamming(), motif, sequence[i:i+length(motif) - 1]) == mismatches) && push!(mtchs, sequence[i:i+length(motif) - 1])
    end
    return mtchs
end

function mismatch_positions(seq1::String, seq2::String) :: Array{Int, 1}
    out::Array{Int, 1} = []
    for i in 1:length(seq1)
        (seq1[i] != seq2[i]) && push!(out, i)
    end
    return out
end