
"""
    Batch_OTA(species, Gene, patterns, Positions; rg=2:18, max_mismatches=4)

Takes as input an array of species, patterns, and positions and creates OTA excel files for every set of pattern and position using every
species indicated. Length of pattern and position arrays must be the same. Gene input is only used in excel file naming. Range of sequence
can be optionally specified as well as the number of mismatches to include up to. Defaults are positions 2-18 of each pattern and 4 mismatches.
"""
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

"""
    Batch_Counts(species, Gene, patterns; rg=2:18)

Takes as input an array of species and patterns and outputs excel file with off target counts from 0 to 4.
Excludes whatever gene of interest is. Default range for patterns is positions 2-18.
"""
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
