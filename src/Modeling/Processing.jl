
function prep_df_hc(df::DataFrame; k::Int=19, alphabet=['A','C','G','U'], colName::String="Sequence") :: DataFrame
	data = DataFrame()
	n = length(alphabet)
	@assert colName in names(df)
	for x in 1:length(df[!, Symbol(colName)])
		(x!=1) && push!(data, [0 for _ in 1:(k*n)])
		for j in 1:n
			for i in 1:k
				(x==1) && (insertcols!(data, i + ((j-1) * k), Symbol("Pos_$(i)_$(alphabet[j])") => 0))
				(df[!, Symbol(colName)][x][i] == alphabet[j]) && (data[!, Symbol("Pos_$(i)_$(alphabet[j])")][x] = 1)
			end
		end
	end
	data
end

function prep_x_hc(df::DataFrame) :: Array{Array{Int64,1},1}
	[dropdims(permutedims(Vector(df[i,:])), dims=1) for i in 1:size(df, 1)]
end

function test_train_split(x, y; split::Float64=0.8)
	@assert length(x) == length(y)
	mnt = rand(length(x)) .< split
	trainX = x[mnt]
	testX = x[.!mnt]
	trainY = y[mnt]
	testY = y[.!mnt]
	return (trainX, trainY, testX, testY)
end

function prep_hc(df::DataFrame; k::Int=19, alphabet=['A', 'C', 'G', 'U'], xName::String="Sequence", yName::String="Activity", split::Float64=0.8)
	@assert xName in names(df)
	@assert yName in names(df)
	df_prepped = prep_df_hc(df, k=k, alphabet=alphabet, colName=xName)
	x = prep_x_hc(df_prepped)
	y = Vector(df[!, Symbol(yName)])
	return test_train_split(x, y, split=split)
end
