


function load_RefSeqDB()
    global REFSEQDB = CSV.read("$(@__DIR__)/RefSeqDB.csv") |> DataFrame
end

load_RefSeqDB()

function list_species()
    for seq in REFSEQDB.Name
        println(seq)
    end
end

function set_species(spec::String = "Human")
    if spec in REFSEQDB.Name
        global SPECIES = spec
        global LINK = REFSEQDB.Link[REFSEQDB.Name .== spec][1]
        global FILEFORMAT = REFSEQDB.FileFormat[REFSEQDB.Name .== spec][1]
        global NUM = REFSEQDB.Num[REFSEQDB.Name .== spec][1]
        global ENCODING = REFSEQDB.Encoding[REFSEQDB.Name .== spec][1]
    else
        println("$(spec) is not an available species.")
        println("To add $(spec) to available species add entry to $(@__DIR__)/RefSeqDB.csv")
        println("To view available species use function list_species")
    end
end

function load_RefSeq()
    global ALLREFSEQ = try 
        collect(values(BSON.load("$PATH/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.bson")))[1];
    catch
        println("Loading $(SPECIES)_mRNA_allRefSeq.bson failed, replace file with save_RefSeq()")
    end
    global GENETRANSCRIPTS = try
        collect(values(BSON.load("$PATH/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.bson")))[1];
    catch
        println("Loading $(SPECIES)_mRNA_GeneTranscripts.bson failed, replace file with save_RefSeq()")
    end
    global TRANSCRIPTGENE = try
        collect(values(BSON.load("$PATH/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.bson")))[1];
    catch
        println("Loading $(SPECIES)_mRNA_TranscriptGene.bson failed, replace file with save_RefSeq()")
    end
end

"""
    download_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Downloads mRNA reference sequence from ftp://ftp.ncbi.nlm.nih.gov/refseq/H\\_sapiens/mRNA\\_Prot/ to the PATH folder.  Defaults to downloading 8 files.
"""
function download_RefSeq(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
    p = Progress(num[end], 0.1, "Updating Reference Sequence ... ", 50)
    ProgressMeter.update!(p, 1; showvalues = [(:File, "$(path)/$(replace(FILEFORMAT, 'X' => '1'))" )])
    for i in num
        download("$(LINK)$(replace(FILEFORMAT, "X" => "$i"))", "$(path)/$(SPECIES)/Download/$(replace(FILEFORMAT, "X" => "$i"))")
        #Check if copy of file already exists and if it does compare hashes with new download
        newfile = "$(path)/$(SPECIES)/$(replace(FILEFORMAT, "X" => "$i"))"
        downfile = "$(path)/$(SPECIES)/Download/$(replace(FILEFORMAT, "X" => "$i"))"
        if isfile(newfile) 
            fd = gzopen(downfile)
            f = gzopen(newfile)
            rd = hash(read(fd, String))
            rf = hash(read(f, String))
            close(fd)
            close(f)
        else
            mv(downfile, newfile, force = true)
            rd = 1
            rf = 1
        end
        if rd != rf
            olddir = replace(split(Libc.strftime(stat(newfile).mtime), " ")[1], "/" => "-")
            !(isdir("$(path)/$(SPECIES)/RefSeq $(olddir)")) && mkdir("$(path)/$(SPECIES)/RefSeq $(olddir)")
            isfile(newfile) && mv(newfile, "$(path)/$(SPECIES)/RefSeq $(olddir)/$(replace(FILEFORMAT, "X" => "$i"))", force=true)
            mv(downfile, newfile, force = true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_allT.bson") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_allT.bson", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_allT.bson", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_df.csv") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_df.csv", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_df.csv", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.bson") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.bson", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_GeneTranscripts.bson", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.bson") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.bson", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_TranscriptGene.bson", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.bson") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.bson", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_allRefSeq.bson", force=true)
        end
        for file in readdir("$(path)/$(SPECIES)/Download/")
            rm("$(path)/$(SPECIES)/Download/$file")
        end
        if i < num[end]
            ProgressMeter.next!(p; showvalues = [(:File, newfile )])
        else
            ProgressMeter.next!(p; showvalues = [(:File, "Finished!" )])
        end
    end
end


"""
    process_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Processes raw gzipped fasta files into DataFrame and saves it as a CSV in the PATH folder
"""
function process_RefSeq(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
    df = DataFrame(Name=String[], ID=String[], Gene=String[], Variant=UInt8[], Sequence=String[], Type=String[])
    d = Dict()
    d[1] = (i -> [split(i, "RNA\n")[1], split(i, " ")[1], length(collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))[end].match[2:end-2] : "None", 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "RNA\n")[1][1:2]])
    d[2] = (i -> [split(i, "RNA\n")[1], split(i, "|")[4], length(collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))[end].match[2:end-2] : "None", 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "|")[4][1:2]])
    for j in num
        f = gzopen("$(path)/$(SPECIES)/$(replace(FILEFORMAT, "X" => "$j"))")
        iter = split(read(f, String), ">")[2:end]
        p = Progress(length(iter), 0.1, "Processing $(replace(FILEFORMAT, "X" => "$j")) ... ")
        for i in iter
            push!(df, d[ENCODING](i))
            ProgressMeter.next!(p)
        end
        close(f)
    end
    CSV.write("$(path)/$(SPECIES)/$(SPECIES)_mRNA_df.csv", df);
    println("Processed data saved to $(path)/$(SPECIES)/$(SPECIES)_mRNA_df.csv")
end


"""
    save_RefSeq(::String=PATH)

Saves relevant data structures from the processes mRNA reference sequence for use in searches.
- TranscriptGene => dictionary of Transcripts to Genes
- GeneTranscripts => dictionary of Genes to Transcripts
- allT => dictionary of Transcript name to base sequence as String
- allRefSeq => dictionary of Transcript name to base sequence as ReferenceSequence
"""
function save_RefSeq(path::String=PATH)
    println("This may take some time")
    println("--Loading processed data from CSV--")
    df = CSV.read("$(path)/$(SPECIES)/$(SPECIES)_mRNA_df.csv") |> DataFrame
    println("--Saving Transcript -> Gene Dictionary--")
    TranscriptGene = Dict(zip(df.ID, df.Gene))
    @save "$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.bson" TranscriptGene
    println("--Saving Gene -> Transcript Dictionary--")
    GeneTranscripts = Dict{String, Array{String, 1}}()
    for (transcript, gene) in TranscriptGene
        if !(haskey(GeneTranscripts, gene))
            GeneTranscripts[gene] = [transcript]
        else
            push!(GeneTranscripts[gene], transcript)
        end
    end
    @save "$(path)/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.bson" GeneTranscripts
    println("--Saving Transcript -> RNA String Dictionary--")
    allT = Dict(zip(df.ID, df.Sequence))
    @save "$(path)/$(SPECIES)/$(SPECIES)_mRNA_allT.bson" allT
    println("--Saving Transcript -> binary RNA Data Dictionary--")
    allRefSeq = Dict{String, ReferenceSequence}()
    for (k, v) in allT
        allRefSeq[k] = encode_refseq(v)
    end
    @save "$(path)/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.bson" allRefSeq
    println("Finished!")
end
