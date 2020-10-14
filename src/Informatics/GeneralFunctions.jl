
"""
    reverse_complement(::String) :: String

Takes the reverse complement of an RNA string
"""
function reverse_complement(pattern::String)::String
    pattern = uppercase(pattern)
    BASES::Array{Char,1} = ['A', 'C', 'G', 'U', 'N']
    @assert sum([x in BASES for x in pattern]) == length(pattern)
    complements::Array{Char,1} = ['U', 'G', 'C', 'A', 'N']
    b_to_c = Dict{Char, Char}(zip(BASES, complements))
    r_c::String = ""
    for base in pattern
        @inbounds r_c = b_to_c[base] * r_c
    end
    return r_c
end


function HammingDist(s1::String, s2::String) :: Int
    total::Int=0
    diff::Int=0
    if length(s1) == length(s2)
        total = 0
        for i in 1:length(s1)
            (s1[i] != s2[i]) && (total += 1)
        end
    elseif length(s1) > length(s2)
        diff = length(s1) - length(s2)
        total = minimum([HammingDist(s1[i-length(s2)+1:i], s2) for i in length(s2):length(s1)]) + diff
    else
        diff = length(s2) - length(s1)
        total = minimum([HammingDist(s1, s2[i-length(s1)+1:i]) for i in length(s1):length(s2)]) + diff
    end
    total
end



function find_gene(gene::String)::String
    #if SPECIES != find_longname(species)
    #    (set_species(species)) ? load_RefSeq() : return ""
    #end
    if gene in keys(GENETRANSCRIPTS)
        return gene
    elseif lowercase(gene) in keys(GENETRANSCRIPTS)
        return lowercase(gene)
    elseif uppercase(gene) in keys(GENETRANSCRIPTS)
        return uppercase(gene)
    elseif uppercase(gene[1]) * lowercase(gene[2:end]) in keys(GENETRANSCRIPTS)
        return uppercase(gene[1]) * lowercase(gene[2:end])
    end
end

function Gen_Kmer(k::Int, Transcript::String) :: DataFrame
    df = DataFrame(Sequence=String[])
    @assert Transcript in keys(ALLREFSEQ)
    seq::String = decode_refseq(ALLREFSEQ[Transcript])
    for i in 1:length(seq) - k + 1
        push!(df, [reverse_complement(seq[i:i+k-1])])
    end
    df
end

function Spec_Score(motif::String, transcriptmotif::String)::Float64
    d = Dict{UnitRange{Int64}, Float64}()
    d[1:7] = 1
    d[8:9] = 1.2
    d[10:10] = 1.25
    d[11:11] = 1.5
    d[12:17] = 1.9
    score::Float64 = 0
    if length(motif) == length(transcriptmotif) && length(motif) == 17
        for (rg, value) in d
            score += (HammingDist(motif[rg], transcriptmotif[rg]) * value)
        end
    else
        score=100
    end
    score
end

function Mismatch_Locations(motif::String, transcriptmotif::String)::Array{Int,1}
    locs::Array{Int,1} = []
    if length(motif) == length(transcriptmotif)
        for (i, (x, y)) in enumerate(zip(motif, transcriptmotif))
            (x != y) && (push!(locs, i))
        end
    end
    locs
end
