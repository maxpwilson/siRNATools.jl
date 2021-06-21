function find_gene_matches(motif::String, gene::String, args::SpecArgs)::Array{MismatchLocus, 1}
    out::Array{MismatchLocus,1} = []
    pattern::String = motif[args.rg]
    pre::Int = args.rg[1] - 1
    suf::Int = length(motif) - args.rg[end]
    (args.anti) && (pattern = reverse_complement(pattern))
    (args.verbose ==true) && (p = Progress(length(ALLREFSEQ), 0.1, "Searching genome... "))
    for name in GENETRANSCRIPTS[gene]
        T = ALLREFSEQ[name]
        k::Int = length(pattern)
        matches::Array{Tuple{Int, Int}} = motif_to_transcript_match(calculate_Peq(pattern), k, T, args.min_mm + 1)
        for (m, loc) in matches
            pre_cor::Int = (loc+pre-length(motif)+1) < 1 ? 1-(loc+pre-length(motif)+1) : 0
            suf_cor::Int = (loc+pre) > T.length ? (loc+pre-T.length) : 0
            locrg::UnitRange{Int64} = loc+pre-length(motif)+1+pre_cor:loc+pre-suf_cor
            trg::UnitRange{Int64} = args.rg[1]+suf_cor:args.rg[end]-pre_cor
            tmotif::String = reverse_complement(decode_refseq_partial(T, locrg))
            mmloc::MismatchLocus = MismatchLocus(; rg=args.rg, transcriptrg = trg, transcript=name, gene=TRANSCRIPTGENE[name], motif=motif, transcriptmotif=tmotif, location=locrg)
            (mmloc.difference <= args.min_mm) && push!(out, mmloc)
        end
        (args.verbose == true) && ProgressMeter.next!(p)
    end
    out
end

function gene_homology_search(motif::String, gene::String; kw...) :: DataFrame
    Args = SpecArgs(; kw...)
    df = DataFrame(Species=String[], Mismatches=Any[], Location=Any[])
    for species in filter(i->i!="Names.csv", readdir("$(PATH)/$(VERSION)/Organisms/"))
        set_species(species)
        load_RefSeq()
        gn = find_gene(gene)
        if gn != ""
            matches = find_gene_matches(motif, gn, Args)
            for match in matches
                push!(df, [species, match.difference, match.transcriptrg])
            end
        end
    end
    df
end
