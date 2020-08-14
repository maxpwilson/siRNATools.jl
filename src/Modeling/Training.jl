
function modelname(model::Chain{T}, suffix::String="")::String where {T}
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

function train_model(model::Chain{T}, opt, loss_start, xyset, epochs::Int=100; accuracy::String="Pearson", penalty::String="None", suffix::String="", alpha::Float64=1e-3, folder::String="", max_improve::Int=10, save_on_start=true, desc=false, max_desc=5) where {T}
    trainX, trainY, testX, testY = xyset
    p = generate_penalty(model, alpha, penalty)
    acc = generate_accuracy(model, accuracy)
    loss(x, y) = loss_start(x, y) + p()
    mname = modelname(model, suffix)
    best_acc = round(acc(testX, testY), digits=4)
    last_improvement = 0
    (length(folder) > 0) && (folder[1] == "/") && (folder = folder[2:end])
    (length(folder) > 0) && (folder[end] == "/") && (folder = folder[1:end-1])
    epoch = 0
    current_acc = best_acc
    train_acc = round(acc(trainX, trainY), digits=4)
    sum_loss_test = mean(loss.(testX, testY))
    sum_loss_train = mean(loss.(trainX, trainY))
    save_on_start && BSON.@save "$(PATH)/$folder/Models/$(mname).bson" model epoch current_acc
    plt = plot([0], [current_acc], ylabel="Accuracy", label = "Test", legend=:topleft, ylims=(-0.1, 1), title="$mname", xlims=(0, epochs))
    plot!(plt, [0], [train_acc], label = "Train")
    plt2 = plot([0], [sum_loss_test], label = "Test", ylabel ="Loss", xlabel="Epochs", legend=false, xlims=(0, epochs))
    plot!(plt2, [0], [sum_loss_train], label = "Train")
    pltc = plot(plt, plt2, layout = (2, 1), size=(2400,1800))
    display(plt)
    @info("Beginning training loop for $mname")
    @info("Epoch 0 accuracy of $(best_acc), $(train_acc)")
    for epoch in 1:epochs
        Flux.train!(loss, params(model), zip(trainX, trainY), opt)
        current_acc = round(acc(testX, testY), digits=4)
        train_acc = round(acc(trainX, trainY), digits=4)
        sum_loss_test = mean(loss.(testX, testY))
        sum_loss_train = mean(loss.(trainX, trainY))
        push!(plt, 1, epoch, current_acc)
        push!(plt, 2, epoch, train_acc)
        push!(plt2, 1, epoch, sum_loss_test)
        push!(plt2, 2, epoch, sum_loss_train)
        png(pltc, "$(PATH)/$folder/Images/$(mname).png")
        @info("Epoch $(epoch) accuracy of $(current_acc), $(train_acc)")
        if current_acc > best_acc
            @info(" -> New best accuracy!")
            BSON.@save "$(PATH)/$folder/Models/$(mname).bson" model epoch current_acc
            best_acc = current_acc
            last_improvement = epoch
        end
        if desc && epoch - last_improvement >= max_desc
            if opt.eta > 1e-6
                opt.eta /= 10.0
                vline!(plt, [epoch], line = (:red, 1), label="η $(opt.eta)")
                vline!(plt2, [epoch], line = (:red, 1))
                @info(" -> Haven't improved in $max_desc epochs, lowering learning rate to $(opt.eta)")
                last_improvement = epoch
            else
                opt.eta = 0.001
                vline!(plt, [epoch], line = (:green, 1), label="η $(opt.eta)")
                vline!(plt2, [epoch], line = (:green, 1))
                @info(" -> Haven't improved in $max_desc epochs, raising learning rate to $(opt.eta)")
                last_improvement = epoch
            end
        end
        display(pltc)
        if max_improve != 0 && epoch - last_improvement >= max_improve
            @warn(" -> Calling this converged")
            break
        end
    end
end
