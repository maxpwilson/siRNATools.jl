module Specificity
using CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter, BSON
using BSON: @save, @load

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

if !("Specificity_Path.txt" in readdir("src"))
    touch("src/Specificity_Path.txt")
end

if !isdir(readline(open("src/Specificity_Path.txt")))
    println("Current specified Path does not exist.  Update using siRNATools.Specificity.Update_Path(PATH). Then re-building siRNATools ")
end
function Update_Path(PATH::String) 
    write("src/Specificity_Path.txt", PATH);
end

const RNA_ALPHABET = RNA_Alphabet(['A', 'C', 'G', 'U', 'N'])
const BASES = RNA_ALPHABET.base_to_bit
const BIT_BASES = RNA_ALPHABET.bit_to_base
const PATH = readline(open("src/Specificity_Path.txt"))
const ALLREFSEQ = try 
    collect(values(BSON.load("$PATH\\Human_mRNA_allRefSeq.bson")))[1]
catch
    println("Loading Human_mRNA_allRefSeq.bson failed, replace file with save_RefSeq()")
end
const GENETRANSCRIPTS = try
    collect(values(BSON.load("$PATH\\Human_mRNA_GeneTranscripts.bson")))[1] 
catch
    println("Loading Human_mRNA_GeneTranscripts.bson failed, replace file with save_RefSeq()")
end
const TRANSCRIPTGENE = try
    collect(values(BSON.load("$PATH\\Human_mRNA_TranscriptGene.bson")))[1]
catch
    println("Loading Human_mRNA_TranscriptGene.bson failed, replace file with save_RefSeq()")
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
    download_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Downloads mRNA reference sequence from ftp://ftp.ncbi.nlm.nih.gov/refseq/H\\_sapiens/mRNA\\_Prot/ to the PATH folder.  Defaults to downloading 7 files.
"""
function download_RefSeq(num::UnitRange{Int64} = 1:7, path::String=PATH)
    p = Progress(num[end], 0.1, "Updating Reference Sequence ... ", 50)
    ProgressMeter.update!(p, 1; showvalues = [(:File, "$(path)\\human.1.rna.fna.gz" )])
    for i in num
        download("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.$(i).rna.fna.gz", "$(path)\\Download\\human.$(i).rna.fna.gz")
        #Check if copy of file already exists and if it does compare hashes with new download
        if isfile("$(path)\\human.$(i).rna.fna.gz") 
            fd = gzopen("$(path)\\Download\\human.$(i).rna.fna.gz")
            f = gzopen("$(path)\\human.$(i).rna.fna.gz")
            rd = hash(read(fd, String))
            rf = hash(read(f, String))
            close(fd)
            close(f)
        else
            mv("$(path)\\Download\\human.$(i).rna.fna.gz", "$(path)\\human.$(i).rna.fna.gz", force = true)
            rd = 1
            rf = 1
        end
        if rd != rf
            olddir = replace(split(Libc.strftime(stat("$(path)\\human.$(i).rna.fna.gz").mtime), " ")[1], "/" => "-")
            !(isdir("$(path)\\RefSeq $(olddir)")) && mkdir("$(path)\\RefSeq $(olddir)")
            isfile("$(path)\\human.$(i).rna.fna.gz") && mv("$(path)\\human.$(i).rna.fna.gz", "$(path)\\RefSeq $(olddir)\\human.$(i).rna.fna.gz", force=true)
            mv("$(path)\\Download\\human.$(i).rna.fna.gz", "$(path)\\human.$(i).rna.fna.gz", force = true)
            isfile("$(path)\\Human_mRNA_allT.bson") && mv("$(path)\\Human_mRNA_allT.bson", "$(path)\\RefSeq $(olddir)\\Human_mRNA_allT.bson", force=true)
            isfile("$(path)\\Human_mRNA_df.csv") && mv("$(path)\\Human_mRNA_df.csv", "$(path)\\RefSeq $(olddir)\\Human_mRNA_df.csv", force=true)
            isfile("$(path)\\Human_mRNA_GeneTranscripts.bson") && mv("$(path)\\Human_mRNA_GeneTranscripts.bson", "$(path)\\RefSeq $(olddir)\\Human_mRNA_GeneTranscripts.bson", force=true)
            isfile("$(path)\\Human_mRNA_TranscriptGene.bson") && mv("$(path)\\Human_mRNA_TranscriptGene.bson", "$(path)\\RefSeq $(olddir)\\Human_mRNA_TranscriptGene.bson", force=true)
            isfile("$(path)\\Human_mRNA_allRefSeq.bson") && mv("$(path)\\Human_mRNA_allRefSeq.bson", "$(path)\\RefSeq $(olddir)\\Human_mRNA_allRefSeq.bson", force=true)
        end
        for file in readdir("$(path)\\Download\\")
            rm("$(path)\\Download\\$file")
        end
        if i < num[end]
            ProgressMeter.next!(p; showvalues = [(:File, "$(path)\\human.$(i+1).rna.fna.gz" )])
        else
            ProgressMeter.next!(p; showvalues = [(:File, "Finished!" )])
        end
    end
end

"""
    process_RefSeq(::UnitRange{Int64}=1:7, ::String=PATH)

Processes raw gzipped fasta files into DataFrame and saves it as a CSV in the PATH folder
"""
function process_RefSeq(num::UnitRange{Int64} = 1:7, path::String=PATH)
    df = DataFrame(Name=String[], ID=String[], Gene=String[], Variant=UInt8[], Sequence=String[], Type=String[])
    for j in num
        f = gzopen("$path\\human.$j.rna.fna.gz")
        iter = split(read(f, String), ">")[2:end]
        p = Progress(length(iter), 0.1, "Processing human.$j.rna.fna.gz ... ")
        for i in iter
            push!(df, [split(i, "RNA\n")[1], split(i, " ")[1], length(collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))[end].match[2:end-2] : "None", 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "RNA\n")[1][1:2]])
            ProgressMeter.next!(p)
        end
        close(f)
    end
    CSV.write("$(path)\\Human_mRNA_df.csv", df);
    println("Processed data saved to $(path)\\Human_mRNA_df.csv")
end

"""
    save_RefSeq(::String=PATH)

Saves relevant data structures from the processes mRNA reference sequence for use in searches.
- TranscriptGene => dictionary of Transcripts to Genes
- GeneTranscripts => dictionary of Genes to Transcripts
- allT => dictionary of Transcript name to base sequence as String
- allRefSeq => dictionary of Transcript name to base sequence as ReferenceSequence
"""
function save_RefSeq(path::String=PATH)
    println("This may take some time")
    println("--Loading processed data from CSV--")
    df = CSV.read("$(path)\\Human_mRNA_df.csv") |> DataFrame
    println("--Saving Transcript -> Gene Dictionary--")
    TranscriptGene = Dict(zip(df.ID, df.Gene))
    @save "$(path)\\Human_mRNA_TranscriptGene.bson" TranscriptGene
    println("--Saving Gene -> Transcript Dictionary--")
    GeneTranscripts = Dict{String, Array{String, 1}}()
    for (transcript, gene) in TranscriptGene
        if !(haskey(GeneTranscripts, gene))
            GeneTranscripts[gene] = [transcript]
        else
            push!(GeneTranscripts[gene], transcript)
        end
    end
    @save "$(path)\\Human_mRNA_GeneTranscripts.bson" GeneTranscripts
    println("--Saving Transcript -> RNA String Dictionary--")
    allT = Dict(zip(df.ID, df.Sequence))
    @save "$(path)\\Human_mRNA_allT.bson" allT
    println("--Saving Transcript -> binary RNA Data Dictionary--")
    allRefSeq = Dict{String, ReferenceSequence}()
    for (k, v) in allT
        allRefSeq[k] = encode_refseq(v)
    end
    @save "$(path)\\Human_mRNA_allRefSeq.bson" allRefSeq
    println("Finished!")
end

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
function find_genome_matches(pattern::String, excluded_gene::String = "",  verbose::Bool = true, minimum_matches = 5) :: Array{Tuple{String, Int64}}
    (length(ALLREFSEQ) == 0) && return []
    out::Array{Tuple{String, Int64}} = []
    (verbose ==true) && (p = Progress(length(ALLREFSEQ), 0.1, "Searching Genome ... "))
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
    df = DataFrame(Pattern=String[], Zero=Int64[], One=Int64[], Two=Int64[], Three=Int64[], Four=Int64[], Score=Float64[])
    for pattern in patterns
        RP = reverse_complement(pattern[rg])
        raw_data = find_genome_matches(RP, excluded_gene, verbose)
        compressed_data = compress_genome_matches(raw_data)
        (mismatchs, spec_score) = final_calc(RP, raw_data, compressed_data)
        push!(df, [pattern, mismatchs[0], mismatchs[1], mismatchs[2], mismatchs[3], mismatchs[4], spec_score])
        (verbose == true) && (println())
    end
    df
end

export Calculate_Specificity, ReferenceSequence

end