using CSV, DataFrames, StatsBase, StringDistances
using BSON: @save, @load

function reverse_complement(pattern::String) :: String
    pattern = uppercase(pattern)
    bases = ['A', 'C', 'G', 'U']
    @assert sum([x in bases for x in pattern]) == length(pattern)
    complements = ['U', 'G', 'C', 'A']
    b_to_c = Dict{String, String}(zip(bases, complements))
    r_c = ""
    for base in pattern
        r_c = b_to_c[base] * r_c
    end
    r_c
end