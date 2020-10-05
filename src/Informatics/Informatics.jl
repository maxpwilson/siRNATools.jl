module Informatics
using CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter, Base.Threads, JuliaDB, MemPool, HTTP, PyCall, FTPClient

include("Path.jl")
include("RNAAlphabet.jl")
include("RefSeq.jl")
include("Maintenance.jl")
include("dbENSEMBL.jl")
include("dbMIRBASE.jl")
include("dbNCBI.jl")
include("GeneralFunctions.jl")
include("SNPs.jl")
include("Specificity.jl")
include("miRNA.jl")
include("BatchProcesses.jl")
include("NameConversion.jl")
include("Transcripts.jl")
include("CreateXlsx.jl")
include("GeneralReport.jl")

export Calculate_Specificity, Deep_Search, getSnps, set_species, list_species, load_RefSeq, Batch_OTA, Batch_Counts

end
