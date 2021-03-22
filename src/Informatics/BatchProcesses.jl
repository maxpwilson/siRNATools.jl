
"""
    Batch_OTA(species, Gene, patterns, Positions; rg=2:18, max_mismatches=4)

Takes as input an array of species, patterns, and positions and creates OTA excel files for every set of pattern and position using every
species indicated. Length of pattern and position arrays must be the same. Gene input is only used in excel file naming. Range of sequence
can be optionally specified as well as the number of mismatches to include up to. Defaults are positions 2-18 of each pattern and 4 mismatches.
"""
function Batch_OTA(Gene::String, Positions::Array{Int,1}; kw...)
    set_species()
    load_RefSeq()
    Patterns = Array{String,1}()
    for p in Positions
        push!(Patterns, Program_Position(Gene, p))
    end
    Batch_OTA(Gene, Patterns, Positions; kw...)
end
#a
function Batch_OTA(Gene::String, Patterns::Array{String,1}, Positions::Array{Int,1}; message::Bool=true, kw...)
    Args = SpecArgs(; kw...)
    @assert length(Patterns) == length(Positions)
    dfs = Dict{Int, Array{DataFrame,1}}()
    for i in 1:length(Patterns)
        dfs[i] = []
    end
    for spec in Args.species
        if set_species(spec)
            load_RefSeq()
            for i in 1:length(Patterns)
                df::DataFrame = Deep_Search(Patterns[i], rg=Args.rg, min_mm=Args.min_mm, anti=Args.anti, verbose=Args.verbose)
                sort!(df, [:MM, :GeneID])
                push!(dfs[i], df)
            end
        end
    end
    for i in 1:length(Patterns)
        dfs[i] = homology_processing(dfs[i], Args.species)
        dfs[i] = expression_processing(dfs[i], Args.species, Args.expression)
        OTA_Excelfile(dfs[i], Args.species, "$(PATH)/Output_Files/OTA$(Args.min_mm)_$(Gene)_Pos$(Positions[i])_rel$(VERSION).xlsx")
        if message
            try
                sendUpdate("Finished OTA$(Args.min_mm)_$(Gene)_Pos$(Positions[i]) for $(join(Args.species, ","))")
            catch
                println("Message Failed")
            end
        end
    end
end

"""
    Batch_Counts(species, Gene, patterns; rg=2:18)

Takes as input an array of species and patterns and outputs excel file with off target counts from 0 to 4.
Excludes whatever gene of interest is. Default range for patterns is positions 2-18.
"""
function Batch_Counts(Gene::String, Patterns::Array{String,1}; kw...)
    Args = SpecArgs(; kw...)
    (Args.excluded_gene == "") && (Args.excluded_gene = Gene)
    dfs::Array{DataFrame,1} = []
    for spec in Args.species
        if set_species(spec)
            load_RefSeq()
            df::DataFrame = Calculate_Specificity(Patterns, excluded_gene=find_gene(Args.excluded_gene), rg=Args.rg, anti=Args.anti)
            push!(dfs, df)
            CSV.write("$(PATH)/Output_Files/$(spec)_Counts_$(Gene).csv", df)
        end
    end
    Counts_Excelfile(dfs, Args.species, "$(PATH)/Output_Files/$(Gene)_Counts.xlsx")
    sendUpdate("Finished $(Gene)_Counts (Gene Length $(length(Patterns))) for $(join(Args.species, ","))")
end
