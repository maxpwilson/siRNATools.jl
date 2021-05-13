function getSNP_CSV(AccID, k=21)
	df = getSnps(AccID, k)
	CSV.write("$(PATH)/Output_Files/SNP_$(AccID).csv", df)
end

function getSNP_Excel(AccID, k=21)
	df = getSnps(AccID, k)
	SNP_ExcelFile(df, "SNP_$(TRANSCRIPTGENE[AccID]).xlsx")
    sendUpdate("Finished SNP_$(TRANSCRIPTGENE[AccID])")
end

function otasnps(df::DataFrame) :: DataFrame
	chr_all = CSV.read("$(PATH)/chrAll.csv", DataFrame)
	df.Pos = [parse(Int,split(i, ":")[1]):parse(Int,split(i, ":")[2]) for i in df[:, "Transcript Location"]]
	genesnps = Dict()
	df.SNP_IDs = ["" for _ in df.Pos]
	df.SNP_Pos = ["" for _ in df.Pos]
	df.SNP_Freq = ["" for _ in df.Pos]
	df.SNP_Change = ["" for _ in df.Pos]
	df.SNP_IDs_in = ["" for _ in df.Pos]
	df.SNP_Pos_in = ["" for _ in df.Pos]
	df.SNP_Freq_in = ["" for _ in df.Pos]
	df.SNP_Change_in = ["" for _ in df.Pos]
	for x in 1:size(df, 1)
		GeneNum = df[x, "Gene ID"]
		GeneName = df[x, "Gene Symbol"]
		AccID = df[x, "Acc#"]
		if !(GeneName in keys(genesnps))
			try
				genesnps[GeneName] = CSV.read("$(PATH)/SNPs/$(SNP_VERSION)/$(GeneName):$(GeneNum).csv", DataFrame)
			catch
				continue
			end
		end
		gene_snps = genesnps[GeneName]
		if !(AccID in chr_all.Transcript)
			continue
		end
		str_rgs = chr_all[chr_all.Transcript.==AccID, :].ChrRange[1]
		rgs = StringToRgs(str_rgs)
		gene_snps_in = gene_snps[in_transcript.(gene_snps.Pos, [rgs for _ in gene_snps.Pos]), :]
		gene_snps_in.TranscriptPos = transcript_position.(gene_snps_in.Pos, [rgs for _ in gene_snps_in.Pos])
		gene_snps_in.Freq = [occursin("1000Genomes:", i) ? parse(Float64, match(r"1000Genomes:[\d\.\,]*\d[\d\.\,]*,(\d[\d\.]*)", i)[1]) : 0 for i in gene_snps_in.Info]
		gene_snps_in.Freq = [i > 0.5 ? 1 - i : i for i in gene_snps_in.Freq]
		gene_snps_in_x = gene_snps_in[[i in df.Pos[x] for i in gene_snps_in.TranscriptPos], :]
		sort!(gene_snps_in_x, [:Pos], rev = !RevRgs(str_rgs))
		gene_snps_in_x.Alt = replace.(gene_snps_in_x.Alt, "T" => "U")
		gene_snps_in_x.Ref = replace.(gene_snps_in_x.Ref, "T" => "U")
		if !RevRgs(str_rgs)
			gene_snps_in_x.Alt = [join(reverse_complement.(split(i, ",")), ",") for i in gene_snps_in_x.Alt]
			gene_snps_in_x.Ref = [join(reverse_complement.(split(i, ",")), ",") for i in gene_snps_in_x.Ref]
		end
		df[x, :].SNP_IDs = join(gene_snps_in_x.ID, ",")
		df[x, :].SNP_Pos = join([21 - (i - df.Pos[x][1]) for i in gene_snps_in_x.TranscriptPos], ",")
		df[x, :].SNP_Freq = join(gene_snps_in_x.Freq, ";")
		df[x, :].SNP_Change = join(["$(gene_snps_in_x[i, :Ref])=>$(gene_snps_in_x[i, :Alt])" for i in 1:size(gene_snps_in_x,1)] , ";")
		mmpos = (!ismissing(df[x, "Mismatch Positions in AS"])) ? parse.(Int, split(df[x, "Mismatch Positions in AS"], ",")) : [0]
		if df[x, :SNP_Pos] == ""
			continue
		end
		for y in 1:length(split(df[x, :SNP_Pos], ","))
			pos = parse.(Int,split(df[x, :SNP_Pos], ","))[y]
			id = split(df[x, :SNP_IDs], ",")[y]
			freq = split(df[x, :SNP_Freq], ";")[y]
			chg = split(df[x, :SNP_Change], ";")[y]
			if pos in mmpos
				df.SNP_IDs_in[x] = df.SNP_IDs_in[x] == "" ? "$(id)" : "$(df.SNP_IDs_in[x]),$(id)"
				df.SNP_Pos_in[x] = df.SNP_Pos_in[x] == "" ? "$(pos)" : "$(df.SNP_Pos_in[x]),$(pos)"
				df.SNP_Freq_in[x] = df.SNP_Freq_in[x] == "" ? "$(freq)" : "$(df.SNP_Freq_in[x]);$(freq)"
				df.SNP_Change_in[x] = df.SNP_Change_in[x] == "" ? "$(chg)" : "$(df.SNP_Change_in[x]);$(chg)"
			end
		end
	end
	df
end

function getSnps(AccID::String, k::Int = 21)
	df = Gen_Kmer(k, AccID)
	df.Pos = [i:i+20 for i in 1:length(df.Sequence)]
	GeneName = ""
	GeneNum = 0
	GeneName = TRANSCRIPTGENE[AccID]
	@assert GeneName != ""
	for j in readdir("$(PATH)/SNPs/$(SNP_VERSION)/")
		if occursin(Regex("^$(GeneName):"), j)
			GeneNum = parse(Int, match(r":([0-9]*)", j)[1])
		end
	end
	@assert GeneNum != 0
	gene_snps = CSV.read("$(PATH)/SNPs/$(SNP_VERSION)/$GeneName:$GeneNum.csv", DataFrame)
	chr_all = CSV.read("$(PATH)/chrAll.csv", DataFrame)
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
function RevRgs_Long(Str::String)::Bool
	occursin("complement", Str)
end


function StringToRgs(Str::String)::Array{Any,1}
    out = []
    for x in split(Str, ",")
        if length(split(x, ":")) > 2
			rg = parse(Int, split(x, ":")[1]):parse(Int, split(x, ":")[2]):parse(Int, split(x, ":")[3])
        	push!(out, rg)
		end
    end
    out
end

function StringToRgs_Long(Str::String)::Array{Any, 1}
	out = []
	c = !occursin("complement", Str)
	stripped_str = replace(Str, r"complement\(|join\(|\)"=>"")
	for rg in split(stripped_str, ",")
		x1 = parse(Int, split(rg, "..")[1])
		x2 = try
			parse(Int, split(rg, "..")[2])
		catch
			x1
		end
		xm = c ? 1 : -1
		rg_out = c ? (x1:xm:x2) : (x2:xm:x1)
		push!(out, rg_out)
	end
	!c && reverse!(out)
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
