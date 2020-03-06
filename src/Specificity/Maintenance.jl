


function load_RefSeqDB()
    try 
        global REFSEQDB = load("$(@__DIR__)/RefSeqDB.jdb") 
    catch
        println("No Reference Sequence Database was found in $(@__DIR__)")
    end
end

load_RefSeqDB()

function list_species()
    for seq in JuliaDB.select(REFSEQDB, :Name)
        println(seq)
    end
end

function set_species(spec::String = "Human")
    if spec in JuliaDB.select(REFSEQDB, :Name)
        global SPECIES = spec
        global LINK = JuliaDB.select(filter(i -> i.Name == spec, REFSEQDB), :Link)[1]
        global FILEFORMAT = JuliaDB.select(filter(i -> i.Name == spec, REFSEQDB), :FileFormat)[1]
        global NUM = JuliaDB.select(filter(i -> i.Name == spec, REFSEQDB), :Num)[1]
        global ENCODING = JuliaDB.select(filter(i -> i.Name == spec, REFSEQDB), :Encoding)[1]
    else
        println("$(spec) is not an available species.")
        println("To add $(spec) to available species add entry to $(@__DIR__)/RefSeqDB.csv")
        println("To view available species use function list_species")
    end
end

function load_RefSeq()
    
    global ALLREFSEQ = try 
        Base.GC.enable(false)
        Dict(load("$PATH/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.jdb"))
        Base.GC.enable(true)
    catch
        println("Loading $(SPECIES)_mRNA_allRefSeq failed, replace file with save_RefSeq()")
    end
    global GENETRANSCRIPTS = try
        Base.GC.enable(false)
        Dict(load("$PATH/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.jdb"))
        Base.GC.enable(true)
    catch
        println("Loading $(SPECIES)_mRNA_GeneTranscripts failed, replace file with save_RefSeq()")
    end
    global TRANSCRIPTGENE = try
        Base.GC.enable(false)
        Dict(load("$PATH/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.jdb"))
        Base.GC.enable(true)
    catch
        println("Loading $(SPECIES)_mRNA_TranscriptGene failed, replace file with save_RefSeq()")
    end
    println("Loaded $SPECIES RefSeq Successfully")
end

function load_RefSeq_old()
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


function full_download_all_RefSeq()
    for name in REFSEQDB.Name
        set_species(name)
        full_download_RefSeq()
    end
end

function full_download_RefSeq()
    println("Fully downloading and processing $(SPECIES) data")
    download_RefSeq()
    println("Download complete, beginning processing")
    process_RefSeq()
    println("Processing complete, beginning save")
    save_RefSeq()
end

"""
    download_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Downloads mRNA reference sequence from ftp://ftp.ncbi.nlm.nih.gov/refseq/H\\_sapiens/mRNA\\_Prot/ to the PATH folder.  Defaults to downloading 8 files.
"""
function download_RefSeq(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
    if !(isdir("$(path)/$(SPECIES)"))
        mkdir("$(path)/$(SPECIES)")
    end
    if !(isdir("$(path)/$(SPECIES)/Download"))
        mkdir("$(path)/$(SPECIES)/Download")
    end
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
            rd = hash(fd)
            rf = hash(f)
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
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_table.jdb") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_table.jdb", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_table.jdb", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.jdb") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.jdb", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_GeneTranscripts.jdb", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.jdb") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.jdb", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_TranscriptGene.jdb", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.jdb") && mv("$(path)/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.jdb", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_allRefSeq.jdb", force=true)
            isfile("$(path)/$(SPECIES)/$(SPECIES)_TranscriptData.jdb") && mv("$(path)/$(SPECIES)/$(SPECIES)_TranscriptData.jdb", "$(path)/$(SPECIES)/RefSeq $(olddir)/$(SPECIES)_mRNA_TranscriptData.jdb", force=true)
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
    tbl = table((Transcript=String[], Range=UnitRange{Int64}[], Type=String[], Gene=String[], GeneID=Int[], Description=String[], RefSeq=ReferenceSequence[]))
    for j in num
        f = gzopen("$(path)/$(SPECIES)/$(replace(FILEFORMAT, "X" => "$j"))")
        p = ProgressUnknown("Processing $(replace(FILEFORMAT, "X" => "$j")) ... ")
        while(!eof(f))
            x = readuntil(f, "LOCUS   ")
            (x == "") && (x = readuntil(f, "LOCUS   "))
            v = match(r"VERSION[\s]*([NX][RM]\_[\d]*\.[\d])", x)[1]
            m = match(r"((\bCDS\b)|(\bncRNA\b)|(\bmisc_RNA\b)|(\brRNA\b))[\s]*(([\<\>\d]*\.\.[\<\>\d]*)|(join\([\<\>\d]*\.\.[\<\>\d]*\,[\<\>\d]*\.\.[\<\>\d]*))", x)
            seq = replace(replace(replace(replace(uppercase(match(r"ORIGIN([\s\S]*)\/\/", x)[1]), " " => ""), "\n" => ""), r"[\d]*" => ""), "T" => "U")
            g = match(r"\/gene=\"(\S*)\"", x)[1]
            d = replace(match(r"\/product=\"([^\"]*)\"", x)[1], "\n                    " => "")
            id = parse(Int, match(r"\"GeneID:([\d]*)\"", x)[1])
            c = m[6]
            t = m[1]
            c = replace(replace(replace(replace(c, ">" => ""), "<" => ""), "," => ".."), "join(" => "")
            cs = split(c, "..")
            c1 = parse(Int, cs[1])
            c2 = parse(Int, cs[end])
            r = c1:c2
            push!(rows(tbl), (Transcript=v, Range=r, Type=t, Gene=g, GeneID=id, Description=d, RefSeq=encode_refseq(seq)))
            ProgressMeter.next!(p)
        end
        ProgressMeter.finish!(p)
        close(f)
    end
    Base.GC.enable(false)
    save(tbl, "$(path)/$(SPECIES)/$(SPECIES)_mRNA_table.jdb")
    Base.GC.enable(true)
    println("Processed data saved to $(path)/$(SPECIES)/$(SPECIES)_mRNA_table")
end


function process_RefSeq_old(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
    println("Warning Deprecated, use process_Refseq")
    df = DataFrame(Name=String[], ID=String[], Gene=String[], Variant=UInt8[], Sequence=String[], Type=String[])
    d = Dict()
    d[1] = (i -> [split(i, "RNA\n")[1], split(i, " ")[1], length(collect(eachmatch(r"(\([^\s:]*[A-Z][^\s:]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([^\s:]*[A-Z][^\s:]*\)\,)", i))[end].match[2:end-2] : collect(eachmatch(r"(\([^\s]*[A-Za-z][^\s]*\)\,)", i))[end].match[2:end-2], 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "RNA\n")[1][1:2]])
    d[2] = (i -> [split(i, "RNA\n")[1], split(i, "|")[4], length(collect(eachmatch(r"(\([^\s:]*[A-Z][^\s:]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([^\s:]*[A-Z][^\s:]*\)\,)", i))[end].match[2:end-2] : collect(eachmatch(r"(\([^\s]*[A-Za-z][^\s]*\)\,)", i))[end].match[2:end-2], 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "|")[4][1:2]])
    d[3] = (i -> [split(i, "RNA\n")[1], split(i, "|")[2], length(collect(eachmatch(r"(\([^\s:]*[A-Z][^\s:]*\)\,)", i))) > 0 ? collect(eachmatch(r"(\([^\s:]*[A-Z][^\s:]*\)\,)", i))[end].match[2:end-2] : collect(eachmatch(r"(\([^\s]*[A-Za-z][^\s]*\)\,)", i))[end].match[2:end-2], 1, replace(replace(split(i, "RNA\n")[2], "\n" => ""), "T" => "U"), split(i, "|")[2][1:2]])
     for j in num
        f = gzopen("$(path)/$(SPECIES)/$(replace(FILEFORMAT, "X" => "$j"))")
        iter = split(read(f, String)[2:end], "\n>")
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
    println("--Loading processed data from JDB--")
    Base.GC.enable(false)
    t = load("$(path)/$(SPECIES)/$(SPECIES)_mRNA_table.jdb") 
    Base.GC.enable(true)
    println("--Saving Transcript -> Gene Dictionary--")
    TranscriptGene = JuliaDB.select(t, (:Transcript, :Gene))
    save(TranscriptGene,"$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.jdb")
    println("--Saving Gene -> Transcript Dictionary--")
    GeneTranscripts = Dict{String, Array{String, 1}}()
    for (transcript, gene) in TranscriptGene
        if !(haskey(GeneTranscripts, gene))
            GeneTranscripts[gene] = [transcript]
        else
            push!(GeneTranscripts[gene], transcript)
        end
    end
    Base.GC.enable(false)
    save(table(GeneTranscripts),"$(path)/$(SPECIES)/$(SPECIES)_mRNA_GeneTranscripts.jdb")
    Base.GC.enable(true)
    println("--Saving Transcript -> binary RNA Data Dictionary--")
    allRefSeq = JuliaDB.select(t, (:Transcript, :RefSeq))
    Base.GC.enable(false)
    save(allRefSeq,"$(path)/$(SPECIES)/$(SPECIES)_mRNA_allRefSeq.jdb")
    Base.GC.enable(true)
    println("--Saving Transcript Data--")
    Base.GC.enable(false)
    save(JuliaDB.select(t, (:Gene, :GeneID, :Description, :Transcript, :Range, :Type)), "$(path)/$(SPECIES)/$(SPECIES)_TranscriptData.jdb")
    Base.GC.enable(true)
    println("Finished!")
end
