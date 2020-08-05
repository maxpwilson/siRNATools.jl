
function Batch_OTA(species,Gene, patterns, Positions; rg::UnitRange{Int64}=2:18, max_mismatches::Int64=5)
    for spec in species
        if set_species(spec)
            load_RefSeq()
            for i in 1:length(patterns)
                df = Deep_Search(patterns[i], rg, max_mismatches)
                CSV.write("$(PATH)/Output_Files/$(spec)_OTA$(max_mismatches-1)_$(Gene)_Pos$(Positions[i]).csv", df)
            end
        end
    end
end

function Batch_Counts(species, Gene, patterns; rg::UnitRange{Int64}=2:18)
    for spec in species
        if set_species(spec)
            load_RefSeq()
            df = Calculate_Specificity(patterns, Gene, rg)
            CSV.write("$(PATH)/Output_Files/$(spec)_Counts_$(Gene).csv", df)
        end
    end
end
