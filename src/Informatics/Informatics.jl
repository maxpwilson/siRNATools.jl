module Informatics
using CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter, Base.Threads, JuliaDB, MemPool, HTTP, PyCall, FTPClient, Telegram
using Parameters: @with_kw

include("Bot.jl")
include("Path.jl")
include("RNAAlphabet.jl")
include("RefSeq.jl")
include("Maintenance.jl")
include("dbENSEMBL.jl")
include("dbMIRBASE.jl")
include("dbNCBI.jl")
include("dbSNP.jl")
include("GeneralFunctions.jl")
include("Structs.jl")
include("SNPs.jl")
include("Search.jl")
include("Program.jl")
include("Specificity.jl")
include("miRNA.jl")
include("OTA_Processing.jl")
include("BatchProcesses.jl")
include("NameConversion.jl")
include("Transcripts.jl")
include("CreateXlsx.jl")
include("GeneralReport.jl")
include("Annotations.jl")
include("HomologySearch.jl")
include("RNAfold.jl")


export Calculate_Specificity, Deep_Search, getSnps, set_species, list_species,
    load_RefSeq, Batch_OTA, Batch_Counts, Gen_Kmer, reverse_complement, download_RefSeq,
    process_RefSeq, Check_NCBI_Version, Check_ENSEMBL_Version, separate_chromosomes,
    annotate_genome, create_chrAll, specialized_annotations, download_ensembl, process_ensembl
    

end
