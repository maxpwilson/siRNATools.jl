function nsnp_processing(dfs :: Array{DataFrame, 1}, species) :: Array{DataFrame, 1}
    dfs_out = Array{DataFrame, 1}()
    for (i, spec) in enumerate(species)
        df = dfs[i]
		missing_transcripts = []
		failed_to_load = []
        if find_longname(spec) == "Homo sapiens"
            chr_all = CSV.read("$(PATH)/chrAll.csv", DataFrame)
            genesnps = Dict()
            df.Neutralizing_ID = ["" for _ in 1:size(df, 1)]
            df.Neutralizing_Pos = ["" for _ in 1:size(df, 1)]
            df.Neutralizing_Freq = ["" for _ in 1:size(df, 1)]
            df.Neutralizing_Change = ["" for _ in 1:size(df, 1)]
			p = Progress(size(df,1), 0.1, "Finding Neutralizing SNPs ...")
            for x in 1:size(df, 1)
                GeneNum = df[x, :GeneID]
                GeneName = df[x, :GeneSymbol]
                AccID = df[x, :Acc]
                if !(GeneName in keys(genesnps))
                    try
                        genesnps[GeneName] = CSV.read("$(PATH)/SNPs/$(SNP_VERSION)/$(GeneName):$(GeneNum).csv", DataFrame)
                    catch
						push!(failed_to_load, GeneName)
                        continue
                    end
                end
                gene_snps = genesnps[GeneName]
                if !(AccID in chr_all.Transcript)
					push!(missing_transcripts, AccID)
                    continue
                end
				try
	                str_rgs = chr_all[chr_all.Transcript .== AccID, :].ChrRange[1]
	                rgs = StringToRgs(str_rgs)
					gene_snps_in = gene_snps[in_transcript.(gene_snps.Pos, [rgs for _ in gene_snps.Pos]), :]
					gene_snps_in.TranscriptPos = transcript_position.(gene_snps_in.Pos, [rgs for _ in gene_snps_in.Pos])
					gene_snps_in.Freq = [occursin("1000Genomes:", i) ? parse(Float64, match(r"1000Genomes:[\d\.\,]*\d[\d\.\,]*,(\d[\d\.]*)", i)[1]) : 0 for i in gene_snps_in.Info]
					gene_snps_in.Freq = [i > 0.5 ? 1 - i : i for i in gene_snps_in.Freq]
					gene_snps_in_strand = gene_snps_in[[i in df.TranscriptLocation[x] for i in gene_snps_in.TranscriptPos], :]
					sort!(gene_snps_in_strand, [:Pos], rev = !RevRgs(str_rgs))
					gene_snps_in_strand.Alt = replace.(gene_snps_in_strand.Alt, "T" => "U")
					gene_snps_in_strand.Ref = replace.(gene_snps_in_strand.Ref, "T" => "U")
					if !RevRgs(str_rgs)
						gene_snps_in_strand.Alt = [join(reverse_complement.(split(i, ",")), ",") for i in gene_snps_in_strand.Alt]
						gene_snps_in_strand.Ref = [join(reverse_complement.(split(i, ",")), ",") for i in gene_snps_in_strand.Ref]
					end
					mmpos = (!ismissing(df[x, :MMPos]) && (df[x,:MMPos]) != "") ? parse.(Int, split(df[x, :MMPos], ",")) : [0]
					gene_snps_in_strand.SNP_Pos = [21 - (i - df.TranscriptLocation[x][1]) for i in gene_snps_in_strand.TranscriptPos]
					gene_snps_neutralizing = gene_snps_in_strand[[(string(gene_snps_in_strand.Alt[i]) == string(uppercase(df.AS[x][gene_snps_in_strand.SNP_Pos[i]]))) &&
					 											(string(gene_snps_in_strand.Ref[i]) == string(uppercase(df.OffTarget[x][gene_snps_in_strand.SNP_Pos[i]]))) &&
																(gene_snps_in_strand.SNP_Pos[i] in mmpos)
																for i in 1:size(gene_snps_in_strand, 1)] , :]

					df[x, :].Neutralizing_ID = join(gene_snps_neutralizing.ID, ",")
					df[x, :].Neutralizing_Pos = join(gene_snps_neutralizing.SNP_Pos, ",")
					df[x, :].Neutralizing_Freq = join(gene_snps_neutralizing.Freq, ";")
					df[x, :].Neutralizing_Change = join(["$(gene_snps_neutralizing[i, :Ref])=>$(gene_snps_neutralizing[i, :Alt])" for i in 1:size(gene_snps_neutralizing,1)] , ";")
				catch
                    print("Failed at position $(x)")
				end
				ProgressMeter.next!(p)
            end
        end
		for x in unique(missing_transcripts)
			println("Transcript $(x) not found in annotations")
		end
		println("Failed to load snps for $(length(unique(failed_to_load))) genes")
	#	for x in unique(failed_to_load)
	#		println("Failed to load snps for $(x)")
	#	end
        push!(dfs_out, df)
    end
    dfs_out
end

function expression_processing(dfs :: Array{DataFrame, 1}, species, expression) :: Array{DataFrame, 1}
    dfs_out = Array{DataFrame, 1}()
    for (i, spec) in enumerate(species)
        df = dfs[i]
        !(spec in keys(expression)) && (expression[spec] = [])
        for tissue in expression[spec]
            try
                df_tissue = CSV.read("$(PATH)/Expression/$(find_longname(spec))/tissue/$(tissue).csv", DataFrame)
                df_tissue = df_tissue[:, vcat(filter(i->occursin(tissue, i), names(df_tissue)), "GeneName")]
                rename!(df_tissue, vcat([i*="_expression" for i in filter(i->occursin(tissue, i), names(df_tissue))] , "GeneName"))
                df = DataFrames.leftjoin(df, df_tissue, on=:GeneSymbol => :GeneName)
                println("$(tissue) expression added for $(spec)")
            catch
                print("Tissue $(tissue) not found")
            end
        end
        push!(dfs_out, df)
    end
    dfs_out
end

function homology_processing(dfs :: Array{DataFrame, 1}, species) :: Array{DataFrame, 1}
    checker = Dict{String, Dict{String, Tuple{Array{Any,1},Array{String,1}}}}()
    for (i, spec) in enumerate(species)
        df = dfs[i]
        d = Dict{String, Tuple{Array{Any,1},Array{String,1}}}()
        for gene in df.GeneSymbol
            d[uppercase(gene)] = (df[df.GeneSymbol .== gene, :].MMPos, replace([k[2:17] for k in df[df.GeneSymbol .== gene, :].OffTarget], r"[A-Z]" => ""))
        end
        checker[spec] = d
    end
    dfs_out = Array{DataFrame, 1}()
    for (i, spec) in enumerate(species)
        df = dfs[i]
        for s in filter(i -> i != spec, species)
            column = [uppercase(i) in keys(checker[s]) ?
            ((length(intersect(checker[s][uppercase(i)][1], checker[spec][uppercase(i)][1])) > 0) ?
            ((length(intersect(checker[s][uppercase(i)][2], checker[spec][uppercase(i)][2])) > 0) ? 3 : 2) : 1) : 0 for i in df.GeneSymbol]
            df[!, "$(s)_homology"] = column
        end
        push!(dfs_out, df)
    end
    return dfs_out
end
