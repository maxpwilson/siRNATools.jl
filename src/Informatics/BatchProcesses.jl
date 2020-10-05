
function Batch_OTA(species,Gene, patterns, Positions; rg::UnitRange{Int64}=2:18, max_mismatches::Int64=4)
    @assert length(patterns) == length(Positions)
    max_mismatches += 1
    dfs = Dict()
    for i in 1:length(patterns)
        dfs[i] = []
    end
    for spec in species
        if set_species(spec)
            load_RefSeq()
            for i in 1:length(patterns)
                df = Deep_Search(patterns[i], rg=rg, max_mismatches=max_mismatches)
                sort!(df, [:MM, :GeneID])
                push!(dfs[i], df)
                #CSV.write("$(PATH)/Output_Files/$(spec)_OTA$(max_mismatches-1)_$(Gene)_Pos$(Positions[i]).csv", df)
            end
        end
    end
    for i in 1:length(patterns)
        OTA_Excelfile(dfs[i], species, "$(PATH)/Output_Files/OTA$(max_mismatches-1)_$(Gene)_Pos$(Positions[i]).xlsx")
    end
end

function Batch_Counts(species, Gene, patterns; rg::UnitRange{Int64}=2:18)
    dfs = []
    for spec in species
        if set_species(spec)
            load_RefSeq()
            df = Calculate_Specificity(patterns, find_gene(Gene, spec), rg)
            push!(dfs, df)
            CSV.write("$(PATH)/Output_Files/$(spec)_Counts_$(Gene).csv", df)
        end
    end
    Counts_Excelfile(dfs, species, "$(PATH)/Output_Files/$(Gene)_Counts.xlsx")
end
