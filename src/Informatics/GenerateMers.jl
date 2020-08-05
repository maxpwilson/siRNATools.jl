
function Gen_Kmer(k::Int, Transcript::String) :: DataFrame
    df = DataFrame(Sequence=String[])
    @assert Transcript in keys(ALLREFSEQ)
    seq = decode_refseq(ALLREFSEQ[Transcript])
    for i in 1:length(seq) - k + 1
        push!(df, [reverse_complement(seq[i:i+k-1])])
    end
    df
end
