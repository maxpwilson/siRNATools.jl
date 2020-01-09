
function load_RefSeq()
    global ALLREFSEQ = try 
        collect(values(BSON.load("$PATH/Human_mRNA_allRefSeq.bson")))[1];
    catch
        println("Loading Human_mRNA_allRefSeq.bson failed, replace file with save_RefSeq()")
    end
    global GENETRANSCRIPTS = try
        collect(values(BSON.load("$PATH/Human_mRNA_GeneTranscripts.bson")))[1];
    catch
        println("Loading Human_mRNA_GeneTranscripts.bson failed, replace file with save_RefSeq()")
    end
    global TRANSCRIPTGENE = try
        collect(values(BSON.load("$PATH/Human_mRNA_TranscriptGene.bson")))[1];
    catch
        println("Loading Human_mRNA_TranscriptGene.bson failed, replace file with save_RefSeq()")
    end
end

"""
    download_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Downloads mRNA reference sequence from ftp://ftp.ncbi.nlm.nih.gov/refseq/H\\_sapiens/mRNA\\_Prot/ to the PATH folder.  Defaults to downloading 8 files.
"""
function download_RefSeq(num::UnitRange{Int64} = 1:8, path::String=PATH)
    p = Progress(num[end], 0.1, "Updating Reference Sequence ... ", 50)
    ProgressMeter.update!(p, 1; showvalues = [(:File, "$(path)/human.1.rna.fna.gz" )])
    for i in num
        download("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.$(i).rna.fna.gz", "$(path)/Download/human.$(i).rna.fna.gz")
        #Check if copy of file already exists and if it does compare hashes with new download
        if isfile("$(path)/human.$(i).rna.fna.gz") 
            fd = gzopen("$(path)/Download/human.$(i).rna.fna.gz")
            f = gzopen("$(path)/human.$(i).rna.fna.gz")
            rd = hash(read(fd, String))
            rf = hash(read(f, String))
            close(fd)
            close(f)
        else
            mv("$(path)/Download/human.$(i).rna.fna.gz", "$(path)/human.$(i).rna.fna.gz", force = true)
            rd = 1
            rf = 1
        end
        if rd != rf
            olddir = replace(split(Libc.strftime(stat("$(path)/human.$(i).rna.fna.gz").mtime), " ")[1], "/" => "-")
            !(isdir("$(path)/RefSeq $(olddir)")) && mkdir("$(path)/RefSeq $(olddir)")
            isfile("$(path)/human.$(i).rna.fna.gz") && mv("$(path)/human.$(i).rna.fna.gz", "$(path)/RefSeq $(olddir)/human.$(i).rna.fna.gz", force=true)
            mv("$(path)/Download/human.$(i).rna.fna.gz", "$(path)/human.$(i).rna.fna.gz", force = true)
            isfile("$(path)/Human_mRNA_allT.bson") && mv("$(path)/Human_mRNA_allT.bson", "$(path)/RefSeq $(olddir)/Human_mRNA_allT.bson", force=true)
            isfile("$(path)/Human_mRNA_df.csv") && mv("$(path)/Human_mRNA_df.csv", "$(path)/RefSeq $(olddir)/Human_mRNA_df.csv", force=true)
            isfile("$(path)/Human_mRNA_GeneTranscripts.bson") && mv("$(path)/Human_mRNA_GeneTranscripts.bson", "$(path)/RefSeq $(olddir)/Human_mRNA_GeneTranscripts.bson", force=true)
            isfile("$(path)/Human_mRNA_TranscriptGene.bson") && mv("$(path)/Human_mRNA_TranscriptGene.bson", "$(path)/RefSeq $(olddir)/Human_mRNA_TranscriptGene.bson", force=true)
            isfile("$(path)/Human_mRNA_allRefSeq.bson") && mv("$(path)/Human_mRNA_allRefSeq.bson", "$(path)/RefSeq $(olddir)/Human_mRNA_allRefSeq.bson", force=true)
        end
        for file in readdir("$(path)/Download/")
            rm("$(path)/Download/$file")
        end
        if i < num[end]
            ProgressMeter.next!(p; showvalues = [(:File, "$(path)/human.$(i+1).rna.fna.gz" )])
        else
            ProgressMeter.next!(p; showvalues = [(:File, "Finished!" )])
        end
    end
end


"""
    process_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Processes raw gzipped fasta files into DataFrame and saves it as a CSV in the PATH folder
"""
function process_RefSeq(num::UnitRange{Int64} = 1:8, path::String=PATH)
    df = DataFrame(Name=String[], ID=String[], Gene=String[], Variant=UInt8[], Sequence=String[], Type=String[])
    for j in num
        f = gzopen("$path/human.$j.rna.fna.gz")
        iter = split(read(f, String), ">")[2:end]
        p = Progress(length(iter), 0.1, "Processing human.$j.rna.fna.gz ... ")
        for i in iter
            push!(df, [split(i, "RNA\n")[1], split(i, " ")[1], length(collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([A-Z][A-Za-z0-9-._/]*\)\,)", i))[end].match[2:end-2] : "None", 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "RNA\n")[1][1:2]])
            ProgressMeter.next!(p)
        end
        close(f)
    end
    CSV.write("$(path)/Human_mRNA_df.csv", df);
    println("Processed data saved to $(path)/Human_mRNA_df.csv")
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
    df = CSV.read("$(path)/Human_mRNA_df.csv") |> DataFrame
    println("--Saving Transcript -> Gene Dictionary--")
    TranscriptGene = Dict(zip(df.ID, df.Gene))
    @save "$(path)/Human_mRNA_TranscriptGene.bson" TranscriptGene
    println("--Saving Gene -> Transcript Dictionary--")
    GeneTranscripts = Dict{String, Array{String, 1}}()
    for (transcript, gene) in TranscriptGene
        if !(haskey(GeneTranscripts, gene))
            GeneTranscripts[gene] = [transcript]
        else
            push!(GeneTranscripts[gene], transcript)
        end
    end
    @save "$(path)/Human_mRNA_GeneTranscripts.bson" GeneTranscripts
    println("--Saving Transcript -> RNA String Dictionary--")
    allT = Dict(zip(df.ID, df.Sequence))
    @save "$(path)/Human_mRNA_allT.bson" allT
    println("--Saving Transcript -> binary RNA Data Dictionary--")
    allRefSeq = Dict{String, ReferenceSequence}()
    for (k, v) in allT
        allRefSeq[k] = encode_refseq(v)
    end
    @save "$(path)/Human_mRNA_allRefSeq.bson" allRefSeq
    println("Finished!")
end
