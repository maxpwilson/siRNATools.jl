function getSnps(AccID::String, k::Int = 21)
	df = Gen_Kmer(k, AccID)
	df.Pos = [i:i+20 for i in 1:length(df.Sequence)]
	GeneName = ""
	GeneNum = 0
	GeneName = TRANSCRIPTGENE[AccID]
	@assert GeneName != ""
	for j in readdir("$(PATH)/SNPs")
		if occursin("$(GeneName):", j)
			GeneNum = parse(Int, match(r":([0-9]*)", j)[1])
		end
	end
	@assert GeneNum != 0
	gene_snps = CSV.read("$(PATH)/SNPs/$GeneName:$GeneNum.csv")
	chr_all = CSV.read("$(PATH)/chrAll.csv")
	str_rgs = chr_all[chr_all.Transcript.==AccID, :].ChrRange[1]
	rgs = StringToRgs(str_rgs)
	gene_snps_in = gene_snps[in_transcript.(gene_snps.Pos, [rgs for _ in gene_snps.Pos]), :]
	gene_snps_in.TranscriptPos = transcript_position.(gene_snps_in.Pos, [rgs for _ in gene_snps_in.Pos])
	gene_snps_in.Freq = [occursin("1000Genomes:", i) ? parse(Float64, match(r"1000Genomes:[\d\.\,]*\d[\d\.\,]*,(\d[\d\.]*)", i)[1]) : 0 for i in gene_snps_in.Info]
	gene_snps_in.Freq = [i > 0.5 ? 1 - i : i for i in gene_snps_in.Freq]
	df.SNP_IDs = ["" for _ in df.Pos]
	df.SNP_Pos = ["" for _ in df.Pos]
	df.SNP_Freq = ["" for _ in df.Pos]
	df.SNP_IDs_wFreq = ["" for _ in df.Pos]
	df.SNP_Pos_wFreq = ["" for _ in df.Pos]
	df.SNP_Freq_wFreq = ["" for _ in df.Pos]
	df.SNP_IDs_hFreq = ["" for _ in df.Pos]
	df.SNP_Pos_hFreq = ["" for _ in df.Pos]
	df.SNP_Freq_hFreq = ["" for _ in df.Pos]
	for x in 1:length(df.Pos)
		gene_snps_in_x = gene_snps_in[[i in df.Pos[x] for i in gene_snps_in.TranscriptPos], :]
		sort!(gene_snps_in_x, [:Pos], rev = RevRgs(str_rgs))
		df[x, :].SNP_IDs = join(gene_snps_in_x.ID, ",")
		df[x, :].SNP_Pos = join([i - df.Pos[x][1] + 1 for i in gene_snps_in_x.TranscriptPos], ",")
		df[x, :].SNP_Freq = join(gene_snps_in_x.Freq, ";")
		df[x, :].SNP_IDs_wFreq = join(gene_snps_in_x[gene_snps_in_x.Freq .> 0 , :].ID, ",")
		df[x, :].SNP_Pos_wFreq = join([i - df.Pos[x][1] + 1 for i in gene_snps_in_x[gene_snps_in_x.Freq .> 0 , :].TranscriptPos], ",")
		df[x, :].SNP_Freq_wFreq = join(gene_snps_in_x[gene_snps_in_x.Freq .> 0 , :].Freq, ";")
		df[x, :].SNP_IDs_hFreq = join(gene_snps_in_x[gene_snps_in_x.Freq .>= 0.01, :].ID, ",")
		df[x, :].SNP_Pos_hFreq = join([i - df.Pos[x][1] + 1 for i in gene_snps_in_x[gene_snps_in_x.Freq .>= 0.01, :].TranscriptPos], ",")
		df[x, :].SNP_Freq_hFreq = join(gene_snps_in_x[gene_snps_in_x.Freq .>= 0.01, :].Freq, ";")
	end

	df
end

function RevRgs(Str::String)::Bool
	if split(split(Str, ",")[1], ":")[2] == "-1"
		return true
	else
		return false
	end
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
