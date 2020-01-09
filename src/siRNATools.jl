module siRNATools

include("ActivityPredictions.jl")
include("ExcelScrape.jl")
include("Specificity/Specificity.jl")
export sht_df, shts_df, get_df_predictions, make_predictions, Calculate_Specificity, ReferenceSequence

end # module
