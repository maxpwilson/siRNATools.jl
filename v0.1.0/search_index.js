var documenterSearchIndex = {"docs":
[{"location":"man/act_guide/#Activity-Guide","page":"Activity Guide","title":"Activity Guide","text":"","category":"section"},{"location":"man/spec_guide/#Specificity-Guide","page":"Specificity Guide","title":"Specificity Guide","text":"","category":"section"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"How to use specificity module of siRNATools","category":"page"},{"location":"man/spec_guide/#Initial-Setup-/-Updates","page":"Specificity Guide","title":"Initial Setup / Updates","text":"","category":"section"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"This only needs to be done once after installing siRNATools","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"After installing siRNATools the default PATH directory must be set up first using the Update_Path function.  The function takes a string as its argument, and each folder must be separated by '\\\\'. The PATH directory is persistent so only needs to be updated once.  Package must be rebuilt after PATH is updated.  ","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> using siRNATools\n\njulia> siRNATools.Specificity.Update_Path(\"C:\\\\My\\\\PATH\\\\folder\")\n\njulia> using Pkg; Pkg.build(\"siRNATools\")","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"Next the reference sequence must be downloaded using the download_RefSeq function.  All files and sub-directories required will be automatically created in the PATH folder.  The current default is set to 7 files, but this has changed in the past.  You can check the files at ftp://ftp.ncbi.nlm.nih.gov/refseq/Hsapiens/mRNAProt/.  ","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> using siRNATools\n\njulia> siRNATools.Specificity.download_RefSeq()","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"After the reference sequence raw files are downloaded they have to be processed into a single readable CSV file using the process_RefSeq function.  Once again all necessary files and sub-folders are created automatically, and the default number of files is 7.","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> siRNATools.Specificity.process_RefSeq()","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"Finally the correct datatypes have to be saved as their own BSON files by save_RefSeq.  This could take a significant amount of time. ","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> siRNATools.Specificity.save_RefSeq()","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"If all of this is done correctly then when rebuilding siRNATools there should no longer be any error messages when using siRNATools.  If the reference sequence needs to be updated repeat previous steps.  If the ncbi has updated the reference sequence the previous steps will save the old files in a folder with the date they were downloaded and create a new set.","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> using Pkg; Pkg.build(\"siRNATools\")\n\njulia> using siRNATools","category":"page"},{"location":"man/spec_guide/#Calculating-Specificity","page":"Specificity Guide","title":"Calculating Specificity","text":"","category":"section"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"First you will need a CSV file containing all of the AntiSense strand sequences you want to calculate specificity for.  They all must be in a single column with the first row being the label for that column.  Sequences must be 5'->3' and contain only A, C, G, and U characters.  ","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"Next you will load the CSV file into julia using the CSV and DataFrames packages.  In this example the dataframe is saved to variable df and the row Antisense contains the relevant information.  ","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> using CSV, DataFrames\n\njulia> df = CSV.read(\"C:\\\\my\\\\csv\\\\path\\\\file.csv\") |> DataFrame\n100×4 DataFrame\n│ Row  │ Sense                 │ Antisense             │\n│      │ String                │ String                │\n├──────┼───────────────────────┼───────────────────────┼\n│ 1    │ AUGGCCCUCCCGACACCCUCG │ CGAGGGUGUCGGGAGGGCCAU │\n.\n.\n.","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"Finally, use the Calculate_Specificity function to calculate the specificity score and gene mismatch information.  Using df.Antisense will calculate for all strands in the antisense column of our dataframe.  Fewer can be picked out using the syntax df.Antisense[i:j] which will calculate for only strands i through j.  The target gene can also be specified and will then be left out of the calculation.  Here we save the resulting dataframe as df_results.  ","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> using siRNATools\n\njulia> df_results = siRNATools.Specificity.Calculate_Specificity(df.Antisense, \"EXAMPLE1\")\n100×7 DataFrame\n│ Row │ Pattern               │ Zero  │ One   │ Two   │ Three │ Four  │ Score   │\n│     │ String                │ Int64 │ Int64 │ Int64 │ Int64 │ Int64 │ Float64 │\n├─────┼───────────────────────┼───────┼───────┼───────┼───────┼───────┼─────────┤\n│ 1   │ CGAGGGUGUCGGGAGGGCCAU │ 0     │ 3     │ 55    │ 452   │ 2894  │ 1       │\n.\n.\n.","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"The data can then be saved as another CSV file using the CSV package and the dataframe df_results from the previous step.","category":"page"},{"location":"man/spec_guide/","page":"Specificity Guide","title":"Specificity Guide","text":"julia> CSV.write(\"C:\\\\my\\\\results\\\\path\\\\data.csv\", df_results)","category":"page"},{"location":"man/Informatics/batch_index/#Batch-Functions","page":"Batch Functions","title":"Batch Functions","text":"","category":"section"},{"location":"man/Informatics/batch_index/","page":"Batch Functions","title":"Batch Functions","text":"siRNATools.Informatics.Batch_OTA\nsiRNATools.Informatics.Batch_Counts\n","category":"page"},{"location":"man/act_index/#Activity-Functions","page":"Activity Functions","title":"Activity Functions","text":"","category":"section"},{"location":"#siRNATools","page":"Home","title":"siRNATools","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Package with tools for siRNA","category":"page"},{"location":"","page":"Home","title":"Home","text":"Package was made using Julia 1.5.1.  ","category":"page"},{"location":"#Informatics","page":"Home","title":"Informatics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"man/Informatics/batch_index.md\"\n]","category":"page"},{"location":"#Modeling","page":"Home","title":"Modeling","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n]","category":"page"},{"location":"man/spec_index/#Specificity-Functions","page":"Specificity Functions","title":"Specificity Functions","text":"","category":"section"},{"location":"man/spec_index/","page":"Specificity Functions","title":"Specificity Functions","text":"","category":"page"},{"location":"man/spec_index/#Reference-Sequence","page":"Specificity Functions","title":"Reference Sequence","text":"","category":"section"},{"location":"man/spec_index/","page":"Specificity Functions","title":"Specificity Functions","text":"siRNATools.Informatics.ReferenceSequence\n","category":"page"},{"location":"man/spec_index/#siRNATools.Informatics.ReferenceSequence","page":"Specificity Functions","title":"siRNATools.Informatics.ReferenceSequence","text":"Structure for reference sequences.  Compresses RNA data into 2 bits of information from the 8 of a normal character string.  Can only use bases A, C, G, U.\n\nA => 00\nC => 01\nG => 10\nU => 11\n\n\n\n\n\n","category":"type"},{"location":"man/install_guide/#Installation-Guide","page":"Installation Guide","title":"Installation Guide","text":"","category":"section"},{"location":"man/install_guide/","page":"Installation Guide","title":"Installation Guide","text":"First download and install Julia from here if not already installed","category":"page"},{"location":"man/install_guide/","page":"Installation Guide","title":"Installation Guide","text":"Next start Julia and type ] into the command line.  The screen should look like below.","category":"page"},{"location":"man/install_guide/","page":"Installation Guide","title":"Installation Guide","text":"(v1.5) Pkg>\n","category":"page"},{"location":"man/install_guide/","page":"Installation Guide","title":"Installation Guide","text":"Then type the below code to install package","category":"page"},{"location":"man/install_guide/","page":"Installation Guide","title":"Installation Guide","text":"(v1.5) Pkg> add https://github.com/maxpwilson/siRNATools.jl.git\n","category":"page"}]
}
