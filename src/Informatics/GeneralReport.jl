function get_transcripts(species::Array{String, 1}, Gene::String) :: DataFrame
    df = DataFrame(Species=String[], Gene=String[], Transcript=String[], Sequence=String[],BasePairs=Int[], EnsemblMatch=String[])
    for x in species
        if set_species(x, verbose=false)
            load_RefSeq(verbose=false)
            g = find_gene(Gene, x)
            if g in keys(GENETRANSCRIPTS)
                for t in GENETRANSCRIPTS[g]
                    s = decode_refseq(ALLREFSEQ[t])
                    bp = length(s)
                    push!(df, [x, g, t, s, bp, ""])
                end
                dfe = filter(i -> i.Gene == g, ENSEMBLDATA) |> DataFrame
                for i in 1:size(dfe)[1]
                    if size(df[df.Sequence.==dfe[i, :Sequence], :])[1] == 1
                        df[df.Sequence.==dfe[i, :Sequence], :EnsemblMatch] = dfe[i, :Accession]
                    else
                        push!(df, [x, g, dfe[i, :Accession], dfe[i, :Sequence], length(dfe[i, :Sequence]), ""])
                    end
                end
            end
        end
    end
    df
end

function General_Report(Species::String, XSpecies, Gene::String, Transcript::String; k::Int=21, mSpecies=[])
    tSpecies = Species in XSpecies ? XSpecies : vcat(Species, XSpecies)
    (length(mSpecies) == 0) && (mSpecies = tSpecies)
    println("Getting Transcripts")
    transcripts = get_transcripts(tSpecies, Gene)
    set_species(Species, verbose=false)
    load_RefSeq(verbose=false)
    @assert Transcript in keys(TRANSCRIPTGENE)
    @assert TRANSCRIPTGENE[Transcript] == find_gene(Gene, Species)
    cdsrg = JuliaDB.select(filter(i -> i.Transcript == Transcript, TRANSCRIPTDATA), :Range)[1]
    df = Gen_Kmer(k, Transcript)
    df.Region = [i in cdsrg ? "CDS" : i < cdsrg[1] ? "5'UTR" : "3'UTR" for i in 1:length(df.Sequence)]
    df.SenseSequence = reverse_complement.(df.Sequence)
    println("Calculating X-Reactivity")
    for i in 1:size(transcripts)[1]
        df[!, Symbol("$(transcripts[i, :Species]):2-21:$(transcripts[i, :Transcript])")] = Int.(occursin.([j[1:20] for j in df.SenseSequence], transcripts[i, :Sequence]))
        df[!, Symbol("$(transcripts[i, :Species]):2-18:$(transcripts[i, :Transcript])")] = Int.(occursin.([j[4:20] for j in df.SenseSequence], transcripts[i, :Sequence]))
        df[!, Symbol("$(transcripts[i, :Species]):<2:$(transcripts[i, :Transcript])")] = mid_xreact.(df.SenseSequence, transcripts[i, :Sequence])
    end
    println("Calculating miRNA homology")
    for x in mSpecies
        df[!, Symbol("$x:miRNA")] = [0 for _ in 1:length(df.Sequence)]
    end
    for i in 1:size(df)[1]
        df_mirna = get_miRNA_matches(df.Sequence[i])
        for x in mSpecies
            if size(df_mirna[df_mirna.Organism .== find_longname(x), :])[1] > 0
                df[i, Symbol("$x:miRNA")] = 1
            end
        end
    end
    println("Analyzing sequences")
    for (key, value) in seq_evaluators()
        df[!, Symbol(key)] = value.(df.Sequence)
    end
    df.CIscore = CI_score.(df.Sequence)

    println("Done")
    GeneralReport_Excelfile(df, transcripts, tSpecies, mSpecies, Transcript, "$(Gene) Informatics.xlsx", Species)
end

function mid_xreact(motif, sequence)
    if occursin(motif[1:20], sequence)
        return 1
    end
    if !occursin(motif[4:20], sequence)
        return 0
    end
    for rg in findall(motif[4:20], sequence)
        if rg[1] > 3
            if (sequence[rg[1]-1] == motif[3] && sequence[rg[1]-2] == motif[2]) || (sequence[rg[1]-3] == motif[1] && sequence[rg[1]-2] == motif[2]) || (sequence[rg[1]-1] == motif[3] && sequence[rg[1]-3] == motif[1])
                return 1
            end
        end
    end
    return 0
end

function seq_evaluators()::Dict
    d = Dict()
    d["GC 2-21"] = x -> count("G", x[2:21]) + count("C", x[2:21])
    d["GC 2-7"] = x -> count("G", x[2:7]) + count("C", x[2:7])
    d["Δ₄GC"] = x -> (count("G", x[18:21]) + count("C", x[18:21])) - (count("G", x[2:4]) + count("C", x[2:4]))
    d["GC 5-9"] = x -> count("G", x[5:9]) + count("C", x[5:9])
    d["ΔₜGC"] = x -> (count("G", x[8:13]) + count("C", x[8:13])) - (count("G", x[2:7]) + count("C", x[2:7]))
    d["GC 16-18"] = x -> count("G", x[16:18]) + count("C", x[16:18])
    d["GC 19-21"] = x -> count("G", x[19:21]) + count("C", x[19:21])
    d["GC first"] = x -> (occursin("G", x[2:21]) ? findfirst("G", x[2:21])[1] + 1 : 0) <= (occursin("C", x[2:21]) ? findfirst("C", x[2:21])[1] : 0) + 1 ? (occursin("G", x[2:21]) ? findfirst("G", x[2:21])[1] : 0) + 1 : (occursin("C", x[2:21]) ? findfirst("C", x[2:21])[1] + 1 : 0)
    d["4xGC 2-20"] = x -> occursin("GGGG", replace(x[2:20], "C" => "G")) ? 1 : 0
    d["6xAU 8-21"] = x -> occursin("UUUUUU", replace(x[8:21], "A" => "U")) ? 1 : 0
    d["3,5,6"] = x -> (x[3] == 'A' || x[3] == 'G') + (x[5] == 'A' || x[5] == 'G') + (x[6] == 'A' || x[6] == 'G')
    d["7,10,14"] = x -> (x[7] == 'A' || x[7] == 'U') + (x[10] == 'A' || x[10] == 'U') + (x[14] == 'A' || x[14] == 'U')
    d
end

function CI_score(x)
    d = seq_evaluators()
    ci = d["3,5,6"](x)
    ci += (0.25 * (((x[8] == 'A') || (x[8] == 'G')) + ((x[9] == 'A') || (x[9]) == 'G')))
    ci += ((x[7] == 'C' && x[6] != 'A') ? -1 : x[7] == 'G' ? 0 : x[7] == 'U' ? 0.5 : 1)
    ci += (x[10] == 'C' ? -1 : (x[10] == 'G' ? 0 : (x[10] == 'U' ? 0.5 : 1)))
    ci += (x[14] == 'G' ? -0.5 : (x[14] == 'C' ? 0 : (x[14] == 'A' ? 0.5 : 1)))
    ci += (x[17] == 'C' ? -0.5 : (x[17] == 'A' ? 0 : (x[17] == 'G' ? 0.25 : 0.5)))
    ci += (x[15] == 'U' ? -0.5 : (x[15] == 'G' ? 0 : (x[15] == 'A' ? 0.25 : 0.5)))
    ci += (d["ΔₜGC"](x) < -1 ? d["ΔₜGC"](x) : 0)
    ci += (d["GC first"](x) > 4 ? ((((x[2] != 'U') && (x[3] != 'U') && (x[4] != 'U')) || (d["GC 2-7"](x) > 1) || (d["GC 2-21"](x)>11)) ? ((d["GC first"](x) - 4) * -0.5) : 0) : 0)
    ci += (d["4xGC 2-20"](x) * -1)
    ci += (d["6xAU 8-21"](x) * -1)
    ci += (d["GC 2-7"](x) > 3 ? (3 - d["GC 2-7"](x)) : 0)
    ci += ((((d["GC 19-21"](x) == 0) && (d["GC 2-21"](x) > 9)) || ((d["GC 19-21"](x) == 1) && (d["GC 2-21"](x) > 10))) ? -1 : 0)
    ci += (d["GC 5-9"](x) > 3 ? 3 - d["GC 5-9"](x) : 0)
    ci
end
