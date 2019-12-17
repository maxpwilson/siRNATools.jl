module siRNATools

include("ActivityPredictions.jl")
include("ExcelScrape.jl")
include("Specificity.jl")
export sht_df, shts_df, get_df_predictions, make_predictions, Calculate_Specificity, load_RefSeq, unload_RefSeq, ReferenceSequence

end # module
