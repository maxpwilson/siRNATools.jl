using DataFrames, CSV, Flux, Query
using BSON: @save, @load

function get_prediction(model, x::Array{Int64,1}, rnd=true::Bool)
    pred = Flux.data(model(x))[1]
    (rnd == true) && (pred = round(pred))
    pred
end

function get_predictions(model, x::Array{Array{Int64,1},1}, rnd=true::Bool)
    a=[]
    for y in x
        push!(a, get_prediction(model, y, rnd))
    end
    a
end

function get_df_predictions(pred::String, model, df::DataFrame, rnd=true::Bool) :: DataFrame
    x = prep_x(df)
    y_hat = get_predictions(model, x, rnd)
    df[!, Symbol("Pred_$pred")] = y_hat
    df
end
function get_df_predictions(pred::String, model_name::String, df::DataFrame, rnd=true::Bool) :: DataFrame
    model = load_model(model_name)
    get_df_predictions(pred, model, df, rnd)
end
function get_df_predictions(pred::String, model, df_file::String, rnd=true::Bool) :: DataFrame
    df = CSV.read(df_file) |> DataFrame
    get_df_predictions(pred, model, df, rnd)
end
function get_df_predictions(pred::String, model_name::String, df_file::String, rnd=true::Bool) :: DataFrame
    df = CSV.read(df_file) |> DataFrame
    model = load_model(model_name)
    get_df_predictions(pred, model, df, rnd)
end

function prep_x(df::DataFrame) :: Array{Array{Int64,1},1}
    data = prep_df(df)
    [dropdims(permutedims(Vector(data[i, :])), dims=1) for i in 1:size(data,1)]
end

function prep_df(df::DataFrame) :: DataFrame
    data = DataFrame()
    b = ['A', 'C', 'G', 'U']
    for x in 1:length(df.Antisense)
        (x!=1) && push!(data, [0 for _ in 1:(19*4)])
        for j in 1:4
            for i in 1:19
                (x==1) && (insertcols!(data, i + ((j-1)*19), Symbol("Pos_$(i)_$(b[j])") => 0))
                (df.Antisense[x][i] == b[j]) && (data[!, Symbol("Pos_$(i)_$(b[j])")][x] = 1)
            end
        end
    end
    data
end

function load_model(filename::String)
    try
        @load filename model
        return model
    catch
        println("Failed to load model $filename")
    end
end
