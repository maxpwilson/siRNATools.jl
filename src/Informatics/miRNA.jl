
function save_miRNA_matches(pat, name, sense::Bool = false)
    df = get_miRNA_matches(pat, sense)
    miRNA_ExcelFile(df, name)
    #CSV.write("$(PATH)/Output_Files/$(name)", df)
end

function get_miRNA_matches(pat, sense::Bool = false)::DataFrame
    df = MIRBASE |> DataFrame
    sense && (pat = reverse_complement(pat))
    index = [i[2:7] == pat[2:7] for i in df.Seq]
    df = sort!(df[index, :], [:CommonName, :Organism])
    replace!(df.CommonName, missing => "Not Found")
    df
end
