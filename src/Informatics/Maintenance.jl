

function list_species()
    if isdir("$(PATH)/$(VERSION)/Organisms")
        for folder in readdir("$(PATH)/$(VERSION)/Organisms")
            (folder[end-2:end] != "csv") && println("$(folder) - $(split(find_vernacular(folder), ",")[1])")
        end
    end
end

function set_species(spec = "Homo sapiens"; verbose::Bool=true)
    if isdir("$(PATH)/$(VERSION)/Organisms")
        if spec in readdir("$(PATH)/$(VERSION)/Organisms")
            (verbose) && println("$spec set")
            global SPECIES = spec
            true
        elseif find_longname(spec) in readdir("$(PATH)/$(VERSION)/Organisms")
            (verbose) && println("$(find_longname(spec)) set")
            global SPECIES = find_longname(spec)
            true
        elseif length(search_vernacular(spec)) + length(search_longnames(spec)) > 1
            (verbose) && println("$(spec) is ambiguous, the following species are available:")
            printed = []
            for x in search_vernacular(spec)
                (verbose) && println("$(find_longname(x)) - $(split(x, ",")[1])")
                push!(printed, find_longname(x))
            end
            for x in search_longnames(spec)
                !(x in printed) && println("$(x) - $(split(find_vernacular(x), ",")[1])")
            end
            false
        elseif length(search_vernacular(spec)) == 1
            (verbose) && println("$(find_longname(search_vernacular(spec)[1])) set")
            global SPECIES = find_longname(search_vernacular(spec)[1])
            true
        elseif length(search_longnames(spec)) == 1
            (verbose) && println("$(search_longnames(spec)[1]) set")
            global SPECIES = search_longnames(spec)[1]
            true
        else
            (verbose) && println("$(spec) is not a valid species")
            false
        end
    else
        (verbose) && println("RefSeq is not downloaded or processed")
    end
end

function load_RefSeq(;verbose::Bool=true)
    if isdir("$(PATH)/$(VERSION)/Organisms/$(SPECIES)")
        Base.GC.enable(false);
        global ALLREFSEQ = try
            Dict(load("$PATH/$VERSION/Organisms/$SPECIES/allRefSeq.jdb"))
        catch
            (verbose) && println("Loading allRefSeq.jdb failed, replace file with save_RefSeq()")
        end
        global GENETRANSCRIPTS = try
            Dict(load("$PATH/$VERSION/Organisms/$SPECIES/GeneTranscripts.jdb"))
        catch
            (verbose) && println("Loading GeneTranscripts.jdb failed, replace file with save_RefSeq()")
        end
        global TRANSCRIPTGENE = try
            Dict(load("$PATH/$VERSION/Organisms/$SPECIES/TranscriptGene.jdb"))
        catch
            (verbose) && println("Loading TranscriptGene.jdb failed, replace file with save_RefSeq()")
        end
        global TRANSCRIPTDATA = try
            load("$PATH/$VERSION/Organisms/$SPECIES/TranscriptData.jdb")
        catch
            (verbose) && println("Loading TranscriptData.jdb failed.")
        end
        global ENSEMBLDATA = try
            load("$PATH/Ensembl/$(ENSEMBL_VERSION)/Organisms/$(SPECIES).jdb")
        catch
            (verbose) && println("Ensembl data not available for $(SPECIES)")
        end
        Base.GC.enable(true);
        (verbose) && println("Loaded $SPECIES RefSeq Successfully")
        true
    else
        (verbose) && println("Species not properly loaded")
        false
    end
end
