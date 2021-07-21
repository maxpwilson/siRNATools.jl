function Program_Position(gene::String, pos::Int; k::Int=21)::String
    df = CSV.read("$(PATH)/ProgramInfo.csv", DataFrame)
    @assert gene in df.Gene
    acc = df[df.Gene .== gene, :Human_Acc][1]
    pat = Gen_Kmer(k, acc).Sequence[pos]
    return pat
end
function Program_AddAcc(gene::String, acc::String)
    df = CSV.read("$(PATH)/ProgramInfo.csv", DataFrame)
#    print(df)
    push!(df, [gene, acc])
    CSV.write("$(PATH)/ProgramInfo.csv", df)
end
