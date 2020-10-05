function gene_name_variants(gene::String)
    out = []
    push!(out, gene)
    (gene != uppercase(gene)) && push!(out, uppercase(gene))
    (gene != lowercase(gene)) && push!(out, lowercase(gene))
    (gene != uppercase(gene[1]) * lowercase(gene[2:end])) && push!(out, uppercase(gene[1]) * lowercase(gene[2:end]))
    out
end
function get_transcripts(gene::String, species=[])
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
