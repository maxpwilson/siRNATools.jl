
function modelname(model::Chain{T}, suffix::String)::String where {T}
    mname = "NN_HL$(length(model) - 1)"
    mname = mname * "_I$(size(model[1].W, 2))"
    for i in model
        if typeof(i) <: Dense
            mname = mname * "_D$(size(i.W)[1])$(String(Symbol(i.σ))[1])"
        elseif typeof(i) <: Dropout
            mname = mname * "_DP$(i.p)"
        elseif typeof(i) <: AlphaDropout
            mname = mname * "_ADP$(i.p)"
        end
    end
    (length(suffix) > 0 ) && (mname = mname * "_$suffix")
    mname
end

function train_model(model::Chain{T}, opt, loss_start, xyset, epochs::Int=100; accuracy::String="Pearson", penalty::String="None", suffix::String="", alpha::Float64=1e-3, folder::String="", max_improve::Int=10) where {T}
    trainX, trainY, testX, testY = xyset
    p = generate_penalty(model, alpha, penalty)
    acc = generate_accuracy(model, accuracy)
    loss(x, y) = loss_start(x, y)
    mname = modelname(model, suffix)
    best_acc = acc(testX, testY)
    last_improvement = 0
    (length(folder) > 0) && (folder[1] == "/") && (folder = folder[2:end])
    (length(folder) > 0) && (folder[end] == "/") && (folder = folder[1:end-1])

    @info("Beginning training loop for $mname")
    @info("Epoch 0 accuracy of $(best_acc)")
    for epoch in 1:epochs
        Flux.train!(loss, params(model), zip(trainX, trainY), opt)
        current_acc = acc(testX, testY)
        @info("Epoch $(epoch) accuracy of $(current_acc)")
        if current_acc > best_acc
            @info(" -> New best accuracy!")
            BSON.@save "$(PATH)/$folder/$(mname).bson" model epoch acc
            best_acc = current_acc
            last_improvement = epoch
        end
        if max_improve != 0 && epoch - last_improvement >= 10
            @warn(" -> Calling this converged")
            break
        end
    end
end
