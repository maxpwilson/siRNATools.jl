

function gen_snn(hl::Int,inp::Int,out::Int,dropout::Float64,npl::Int)
    layers = []
    push!(layers, Dense(inp, npl, selu))
    push!(layers, AlphaDropout(dropout))
    for i in 1:hl-1
        push!(layers, Dense(npl, npl, selu))
        push!(layers, AlphaDropout(dropout))
    end
    push!(layers, Dense(npl, out, selu))
    Chain(layers...)
end
function gen_snn(hl::Int,inp::Int,out::Int,dropout::Float64,nodestart::Int,nodestop::Int)
    layers = []
    push!(layers, Dense(inp, nodestart, selu))
    push!(layers, AlphaDropout(dropout))
    nodediff = Int(round((nodestart - nodestop) / (hl - 1)))
    currentnode = nodestart
    for i in 1:hl-1
        push!(layers, Dense(currentnode, i == hl-1 ? nodestop : currentnode  - nodediff, selu))
        push!(layers, AlphaDropout(dropout))
        currentnode = currentnode - nodediff
    end
    push!(layers, Dense(nodestop, out, selu))
    Chain(layers...)
end
