function gene_name_variants(gene::String)
    out = []
    push!(out, gene)
    (gene != uppercase(gene)) && push!(out, uppercase(gene))
    (gene != lowercase(gene)) && push!(out, lowercase(gene))
    (gene != uppercase(gene[1]) * lowercase(gene[2:end])) && push!(out, uppercase(gene[1]) * lowercase(gene[2:end]))
    out
end
function get_transcripts_ensembl(gene::String)
    df = DataFrame(ENSEMBLDATA)
    df_out = df[df.Gene .== find_gene(gene), :]
end
function save_transcripts_ensembl(gene::String, species::Array{String,1}=["Human"])
    df = DataFrame(Species=String[], Accession=String[], Gene=String[], Description=String[], Sequence=String[])
    for spec in species
        if set_species(spec)
            load_RefSeq()
            tdf = DataFrame(ENSEMBLDATA)
            tdf[!, :Species] .= spec
            df = vcat(df, tdf[tdf.Gene .== find_gene(gene), :])
        end
    end
    CSV.write("$(PATH)/Output_Files/ENSEMBL_$(gene)_$(join(species, "")).csv", df)
end
function get_transcripts_ncbi(gene::String, species=[])
    if length(species) == 0
        species = readdir("$(PATH)/$(VERSION)/Organisms")
    end
    df = DataFrame(Organism=String[], Transcript=String[])
    for x in species
        if set_species(x)
            load_RefSeq()
            for y in gene_name_variants(gene)
                if y in keys(GENETRANSCRIPTS)
                    for z in GENETRANSCRIPTS[y]
                        push!(df, [x, z])
                    end
                end
            end
        end
    end
    df
end
