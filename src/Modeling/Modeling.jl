module Modeling
using DataFrames, CSV, Flux, Query, Statistics, LinearAlgebra, StatsBase, BSON, Plots

include("Path.jl")
include("ActivityPredictions.jl")
include("Processing.jl")
include("LossAcc.jl")
include("Training.jl")
include("ModelGen.jl")

export gen_snn, train_model, modelname

end
