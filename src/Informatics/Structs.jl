@with_kw mutable struct MismatchLocus
    transcript::String = ""
    gene::String = ""
    motif::String = "A"
    rg::UnitRange{Int64} = 1:1
    transcriptrg::UnitRange{Int64} = rg
    transcriptmotif::String = "A"
    difference::Int = HammingDist(motif[rg], transcriptmotif[transcriptrg])
    location::UnitRange{Int64} = 1:1
    mismatches::Array{Int,1} = Mismatch_Locations(motif[rg], transcriptmotif[transcriptrg]) .+ (rg[1]-1)
    score::Float64 = Spec_Score(motif[rg], transcriptmotif[transcriptrg])
end

@with_kw mutable struct MismatchGene
    gene::String = ""
    loci::Array{MismatchLocus, 1} = []
    difference::Int = (length(loci) > 0) ? minimum([i.difference for i in loci]) : 0
    score::Float64 = (length(loci) > 0) ? minimum([i.score for i in loci]) : 0
end

@with_kw mutable struct SpecArgs
    verbose::Bool = true
    min_mm::Int = 4
    rg::UnitRange{Int64} = 2:18
    excluded_gene::String = ""
    anti::Bool = true
    expression::Dict{String, Array{String,1}} = Dict{String, Array{String,1}}()
    snps::Bool = false
    dense::Bool = false
    species::Array{String, 1} = ["Human", "Cyno", "Rat", "Mouse", "Rhesus"]
end

struct secondary_structure
    pat::String
    folding_sites::String
    free_energy::Float64
    first_bond::Int
    num_bonds::Int
    bonds::Vector{String}
    function secondary_structure(pat)
        folding_output = rnafold(pat)
        new(pat, folding_output[1], folding_output[2], rnafold_firstbond(folding_output[1]), rnafold_numbonds(folding_output[1]), rnafold_bonds(folding_output[1]))
    end
end

function copy(s::SpecArgs)::SpecArgs
    SpecArgs(verbose=s.verbose,min_mm=s.min_mm,rg=s.rg,excluded_gene=s.excluded_gene,anti=s.anti,species=s.species,snps=s.snps,expression=s.expression)
end
