

function find_gene(gene, species)
    if SPECIES != find_longname(species)
        (set_species(species)) ? load_RefSeq() : return ""
    end
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
    seq = decode_refseq(ALLREFSEQ[Transcript])
    for i in 1:length(seq) - k + 1
        push!(df, [reverse_complement(seq[i:i+k-1])])
    end
    df
end
