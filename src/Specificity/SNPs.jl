function getSnps(df::DataFrame, GeneName::String, GeneNum::Int, AccID::String)
	df.Pos = [i:i+20 for i in 1:length(df.Sequence)]
	gene_snps = CSV.read("$(PATH)/SNPs/$GeneName:$GeneNum.csv")
	chr_all = CSV.read("$(PATH)/chrAll.csv")
	str_rgs = chr_all[chr_all.Transcript.==AccID, :].ChrRange[1]
	rgs = StringToRgs(str_rgs)
	gene_snps_in = gene_snps[in_transcript.(gene_snps.Pos, [rgs for _ in gene_snps.Pos]), :]
	gene_snps_in.TranscriptPos = transcript_position.(gene_snps_in.Pos, [rgs for _ in gene_snps_in.Pos])
	df.SNP_IDs = ["" for _ in df.Pos]
	df.SNP_Pos = ["" for _ in df.Pos]
	for x in 1:length(df.Pos)
		gene_snps_in_x = gene_snps_in[[i in df.Pos[x] for i in gene_snps_in.TranscriptPos], :]
		df[x, :].SNP_IDs = join(gene_snps_in_x.ID, ",")
		df[x, :].SNP_Pos = join([i - df.Pos[x][1] + 1 for i in gene_snps_in_x.TranscriptPos], ",")
	end
	df
end

function StringToRgs(Str::String)::Array{Any,1}
    out = []
    for x in split(Str, ",")
        rg = parse(Int, split(x, ":")[1]):parse(Int, split(x, ":")[2]):parse(Int, split(x, ":")[3])
        push!(out, rg)
    end
    out
end

function transcript_position(pos::Int, rg)::Int
	total = 0
	for x in rg
		if pos in x
			(x[1] > pos) ? (total += x[1] - pos + 1) : (total += pos - x[1] + 1)
			break
		else
			total += length(x)
		end
	end
	total
end


function in_transcript(pos::Int, rg)::Bool
	for x in rg
		(pos in x) && (return true)
	end
	false
end
