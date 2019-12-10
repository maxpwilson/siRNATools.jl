using CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter
using BSON: @save, @load

PATH = "C:\\Users\\mwilson\\Notebooks\\Specificity\\"
ALLT = Dict{String, String}()
GENETRANSCRIPTS = Dict{String, Array{String, 1}}()
TRANSCRIPTGENE = Dict{String, String}()

function download_RefSeq(num::UnitRange{Int64} = 1:8, path::String=PATH)
    p = Progress(num[end], 0.1, "Updating Reference Sequence ... ", 50)
    ProgressMeter.update!(p, 1; showvalues = [(:File, "$(path)human.1.rna.fna.gz" )])
    for i in num
        download("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.$(i).rna.fna.gz", "$(path)Download//human.$(i).rna.fna.gz")
        #Check if copy of file already exists and if it does compare hashes with new download
        if isfile("$(path)human.$(i).rna.fna.gz") 
            fd = gzopen("$(path)Download//human.$(i).rna.fna.gz")
            f = gzopen("$(path)human.$(i).rna.fna.gz")
            rd = hash(read(fd, String))
            rf = hash(read(f, String))
            close(fd)
            close(f)
        else
            mv("$(path)Download//human.$(i).rna.fna.gz", "$(path)human.$(i).rna.fna.gz", force = true)
            rd = 1
            rf = 1
        end
        if rd != rf
            olddir = replace(split(Libc.strftime(stat("$(path)human.$(i).rna.fna.gz").mtime), " ")[1], "/" => "-")
            !(isdir("$(path)RefSeq $(olddir)")) && mkdir("$(path)RefSeq $(olddir)")
            isfile("$(path)human.$(i).rna.fna.gz") && mv("$(path)human.$(i).rna.fna.gz", "$(path)RefSeq $(olddir)//human.$(i).rna.fna.gz", force=true)
            mv("$(path)Download//human.$(i).rna.fna.gz", "$(path)human.$(i).rna.fna.gz", force = true)
            isfile("$(path)Human_mRNA_allT.bson") && mv("$(path)Human_mRNA_allT.bson", "$(path)RefSeq $(olddir)//Human_mRNA_allT.bson", force=true)
            isfile("$(path)Human_mRNA_df.csv") && mv("$(path)Human_mRNA_df.csv", "$(path)RefSeq $(olddir)//Human_mRNA_df.csv", force=true)
            isfile("$(path)Human_mRNA_GeneTranscripts.bson") && mv("$(path)Human_mRNA_GeneTranscripts.bson", "$(path)RefSeq $(olddir)//Human_mRNA_GeneTranscripts.bson", force=true)
            isfile("$(path)Human_mRNA_TranscriptGene.bson") && mv("$(path)Human_mRNA_TranscriptGene.bson", "$(path)RefSeq $(olddir)//Human_mRNA_TranscriptGene.bson", force=true)
        end
        for file in readdir("$(path)Download//")
            rm("$(path)Download//$file")
        end
        if i < num[end]
            ProgressMeter.next!(p; showvalues = [(:File, "$(path)human.$(i+1).rna.fna.gz" )])
        else
            ProgressMeter.next!(p; showvalues = [(:File, "Finished!" )])
        end
    end
end

function process_RefSeq(num::UnitRange{Int64} = 1:8, path::String=PATH)
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
    CSV.write("$(path)Human_mRNA_df.csv", df)
end

function save_RefSeq(path::String=PATH)
    df = CSV.read("$(path)Human_mRNA_df.csv") |> DataFrame
    TranscriptGene = Dict(zip(df.ID, df.Gene))
    @save "$(path)Human_mRNA_TranscriptGene.bson" TranscriptGene
    GeneTranscripts = Dict{String, Array{String, 1}}()
    for (transcript, gene) in TranscriptGene
        if !(haskey(GeneTranscripts, gene))
            GeneTranscripts[gene] = [transcript]
        else
            push!(GeneTranscripts[gene], transcript)
        end
    end
    @save "$(path)Human_mRNA_GeneTranscripts.bson" GeneTranscripts
    allT = Dict(zip(df.ID, df.Sequence))
    @save "$(path)Human_mRNA_allT.bson" allT
end

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

function find_genome_matches(pattern::String, excluded_gene::String = "",  verbose::Bool = true, minimum_matches = 5) :: Array{Tuple{String, Int64}}
    (length(ALLT) == 0) && return []
    out::Array{Tuple{String, Int64}} = []
    (verbose ==true) && (p = Progress(length(ALLT), 0.1, "Searching Genome ... "))
    for (name, T) in ALLT
        (excluded_gene != "") && ((name in GENETRANSCRIPTS[excluded_gene]) && continue)
        match::Int64 = motif_to_transcript_match(pattern, T)
        (match < minimum_matches) && push!(out, (name, match))
        (verbose == true) && ProgressMeter.next!(p)
    end
    out
end

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

function mismatch_counts(compressed_data::Dict{String, Array{Int64, 1}}) :: Dict{Int64, Int64}
    out = Dict{Int64, Int64}(zip([0,1,2,3,4], [0,0,0,0,0]))
    for (gene, matches) in compressed_data
        out[minimum(matches)] += 1
    end
    return out
end

function specificity_score(pattern::String, raw_data::Array{Tuple{String, Int64}}) :: Float64
    (length(ALLT) == 0) && return -1
    min_match::Int64 = 5
    for (a, b) in raw_data
        (b < min_match) && (min_match = b)
    end
    min_score::Float64 = 5
    for (name, match) in raw_data
        if minimum(match) == min_match
            match_patterns::Array{String, 1} = find_match_sequences(pattern, ALLT[name], min_match)
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
    (min_score == 2.9) && (min_score = 3)
    return min_score
end

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
        mismatchs = mismatch_counts(compressed_data)
        spec_score = specificity_score(RP, raw_data)
        push!(df, [pattern, mismatchs[0], mismatchs[1], mismatchs[2], mismatchs[3], mismatchs[4], spec_score])
        (verbose == true) && (println())
    end
    df
end

