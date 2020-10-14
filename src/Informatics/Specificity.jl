"""

"""
function find_genome_matches(motif::String, args::SpecArgs)::Array{MismatchLocus,1}
    out::Array{MismatchLocus,1} = []
    pattern::String = motif[args.rg]
    pre::Int = args.rg[1] - 1
    suf::Int = length(motif) - args.rg[end]
    (args.anti) && (pattern = reverse_complement(pattern))
    (args.verbose ==true) && (p = Progress(length(ALLREFSEQ), 0.1, "Searching genome... "))
    for (name, T) in ALLREFSEQ
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
"""

"""
function compress_genome_matches(raw_data::Array{MismatchLocus,1}, args::SpecArgs)::Array{MismatchGene,1}
    uniquegenes::Array{String,1} = unique([i.gene for i in raw_data])
    filter!(i->i != args.excluded_gene, uniquegenes)
    out_size::Int = length(uniquegenes)
    out = Array{MismatchGene,1}(undef, out_size)
    i::Int = 1
    for gn in uniquegenes
        out[i] = MismatchGene(gene=gn, loci=filter(i->i.gene==gn, raw_data))
        i+=1
    end
    out
end

function Deep_Search(pattern::String; kw::Base.Iterators.Pairs...)::DataFrame
    Args::SpecArgs = SpecArgs(; kw...)
    GeneDescription = Dict{String, String}(zip(JuliaDB.select(TRANSCRIPTDATA, :Gene), JuliaDB.select(TRANSCRIPTDATA, :Description)))
    GeneID = Dict{String,Int64}(zip(JuliaDB.select(TRANSCRIPTDATA, :Gene), JuliaDB.select(TRANSCRIPTDATA, :GeneID)))
    TranscriptRange = Dict{String,UnitRange{Int64}}(zip(JuliaDB.select(TRANSCRIPTDATA, :Transcript), JuliaDB.select(TRANSCRIPTDATA, :Range)))
    TranscriptType = Dict{String,String}(zip(JuliaDB.select(TRANSCRIPTDATA, :Transcript), JuliaDB.select(TRANSCRIPTDATA, :Type)))

    df::DataFrame = DataFrame(Acc=String[], GeneID=Any[], GeneSymbol=String[], Description=String[], Region=Any[], MM=Int[], AS=String[], OffTarget=String[], MMPos=Any[], TranscriptLocation=Any[])

    raw_data::Array{MismatchLocus,1} = find_genome_matches(pattern, Args)
    for loc in raw_data
        gid::Any = (loc.gene in keys(GeneID)) ? GeneID[loc.gene] : "na"
        gd::String = (loc.gene in keys(GeneDescription)) ? GeneDescription[loc.gene] : "na"
        tmotif = loc.transcriptmotif
        for k in loc.mismatches
            tmotif = tmotif[1:k-1] * lowercase(tmotif[k]) * tmotif[k+1:end]
        end
        region::Array{String,1} = []
        tstart::Int = (loc.transcript in keys(TranscriptRange)) ? TranscriptRange[loc.transcript][1] : 1
        tstop::Int = (loc.transcript in keys(TranscriptRange)) ? TranscriptRange[loc.transcript][end] : loc.location[end]
        (tstart < loc.location[end] && tstop > loc.location[1]) && (push!(region, (loc.transcript in keys(TranscriptType)) ? TranscriptType[loc.transcript] : "na"))
        (tstart > loc.location[1]) && (push!(region, "5' UTR"))
        (tstop < loc.location[end]) && (push!(region, "3' UTR"))

        push!(df, [loc.transcript, gid, loc.gene, gd, join(region, ","), loc.difference, loc.motif, tmotif, join(string.(loc.mismatches), ","), loc.location])
    end
    df
end

function Calculate_Specificity(patterns; kw...)::DataFrame
    Calculate_Specificity(string.(patterns); kw...)
end
function Calculate_Specificity(patterns::Array{String, 1}; kw...)::DataFrame
    Args::SpecArgs = SpecArgs(; kw...); OutArgs::SpecArgs = copy(Args)
    OutArgs.verbose = false
    nums = Array{String,1}(["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten"])
    df::DataFrame = DataFrame(:ID=>Int[], :Pattern=>String[], [Symbol(i)=>Int64[] for i in nums[1:Args.min_mm+1]]..., :Score=>Float64[])
    counter = Atomic{Int}(0)
    e_gn::String = find_gene(Args.excluded_gene)
    (Args.verbose == true) && ((p = Progress(length(patterns), 0.2, "Calculating Specificity ... ")))
    (Args.verbose == true) && (update!(p, 0))
    @threads for i in 1:length(patterns)
        pattern = patterns[i]
        atomic_add!(counter, 1)
        raw_data::Array{MismatchLocus,1} = find_genome_matches(pattern, OutArgs)
        compressed_data::Array{MismatchGene,1} = compress_genome_matches(raw_data, OutArgs)
        mismatches::Array{Int,1} = []
        for j in 0:Args.min_mm
            push!(mismatches, length(filter(i->i.difference .== j, compressed_data)))
        end
        score::Float64 = minimum([i.score for i in compressed_data])
        push!(df, [i, pattern, mismatches..., score])
        (Args.verbose == true) && ProgressMeter.next!(p)
    end
    sort!(df, :ID)
    df
end

function Calculate_Specificity_Loop()

end
