module Informatics
using CSV, DataFrames, StatsBase, StringDistances, GZip, ProgressMeter, Base.Threads, JuliaDB, MemPool, HTTP

include("Path.jl")
include("RNAAlphabet.jl")
include("RefSeq.jl")
include("Maintenance.jl")
include("GenerateMers.jl")
include("SNPs.jl")
include("Specificity.jl")
include("BatchProcesses.jl")

export Calculate_Specificity, Deep_Search, getSnps, set_species, list_species, load_RefSeq, Batch_OTA, Batch_Counts

end
