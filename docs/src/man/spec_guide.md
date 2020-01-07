# Specificity Guide

How to use specificity module of siRNATools

## Initial Setup / Updates

This only needs to be done once after installing siRNATools

After installing siRNATools the default PATH directory must be set up first using the [`Update_Path`](@ref) function.  The function takes a string as its argument, and each folder must be separated by '\\\\'. The PATH directory is persistent so only needs to be updated once.  Package must be rebuilt after PATH is updated.  
```julia
julia> using siRNATools

julia> siRNATools.Specificity.Update_Path("C:\\My\\PATH\\folder")

julia> using Pkg; Pkg.build("siRNATools")
```

Next the reference sequence must be downloaded using the [`download_RefSeq`](@ref) function.  All files and sub-directories required will be automatically created in the PATH folder.  The current default is set to 7 files, but this has changed in the past.  You can check the files at ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/.  
```julia
julia> using siRNATools

julia> siRNATools.Specificity.download_RefSeq()
```

After the reference sequence raw files are downloaded they have to be processed into a single readable CSV file using the [`process_RefSeq`](@ref) function.  Once again all necessary files and sub-folders are created automatically, and the default number of files is 7.
```julia
julia> siRNATools.Specificity.process_RefSeq()
```

Finally the correct datatypes have to be saved as their own BSON files by [`save_RefSeq`](@ref).  This could take a significant amount of time. 
```julia
julia> siRNATools.Specificity.save_RefSeq()
```

If all of this is done correctly then when rebuilding siRNATools there should no longer be any error messages when using siRNATools.  If the reference sequence needs to be updated repeat previous steps.  If the ncbi has updated the reference sequence the previous steps will save the old files in a folder with the date they were downloaded and create a new set.
```julia
julia> using Pkg; Pkg.build("siRNATools")

julia> using siRNATools
```

## Calculating Specificity

First you will need a CSV file containing all of the AntiSense strand sequences you want to calculate specificity for.  They all must be in a single column with the first row being the label for that column.  Sequences must be 5'->3' and contain only A, C, G, and U characters.  

Next you will load the CSV file into julia using the CSV and DataFrames packages.  In this example the dataframe is saved to variable df and the row Antisense contains the relevant information.  
```julia
julia> using CSV, DataFrames

julia> df = CSV.read("C:\\my\\csv\\path\\file.csv") |> DataFrame
100×4 DataFrame
│ Row  │ Sense                 │ Antisense             │
│      │ String                │ String                │
├──────┼───────────────────────┼───────────────────────┼
│ 1    │ AUGGCCCUCCCGACACCCUCG │ CGAGGGUGUCGGGAGGGCCAU │
.
.
.
```

Finally use the [`Calculate_Specificity`](@ref) function to calculate the specificity score and gene mismatch information.  Using df.Antisense will calculate for all strands in the antisense column of our dataframe.  Fewer can be picked out using the syntax df.Antisense[i:j] which will calculate for only strands i through j.  The target gene can also be specified and will then be left out of the calculation.  Here we save the resulting dataframe as df_results.  
```julia
julia> df_results = siRNATools.Specificity.Calculate_Specificity(df.Antisense, "EXAMPLE1")
100×7 DataFrame
│ Row │ Pattern               │ Zero  │ One   │ Two   │ Three │ Four  │ Score   │
│     │ String                │ Int64 │ Int64 │ Int64 │ Int64 │ Int64 │ Float64 │
├─────┼───────────────────────┼───────┼───────┼───────┼───────┼───────┼─────────┤
│ 1   │ CGAGGGUGUCGGGAGGGCCAU │ 0     │ 3     │ 55    │ 452   │ 2894  │ 1       │
.
.
.
```

The data can then be saved as another CSV file using the CSV package and the dataframe df_results from the previous step.
```julia
julia> CSV.write("C:\\my\\results\\path\\data.csv", df_results)
```