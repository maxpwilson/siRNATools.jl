

function list_species()
    if isdir("$(PATH)/$(VERSION)/Organisms")
        for folder in readdir("$(PATH)/$(VERSION)/Organisms")
            println(folder)
        end
    end
end

function set_species(spec::String = "Homo sapiens")
    if isdir("$(PATH)/$(VERSION)/Organisms")
        if spec in readdir("$(PATH)/$(VERSION)/Organisms")
            global SPECIES = spec
        else
            println("$(spec) is not avalid species")
        end
    else
        println("RefSeq is not downloaded or processed")
    end
end

function load_RefSeq()
    if isdir("$(PATH)/$(VERSION)/Organisms/$(SPECIES)")
        Base.GC.enable(false);
        global ALLREFSEQ = try
            Dict(load("$PATH/$VERSION/Organisms/$SPECIES/allRefSeq.jdb"))
        catch
            println("Loading allRefSeq.jdb failed, replace file with save_RefSeq()")
        end
        global GENETRANSCRIPTS = try
            Dict(load("$PATH/$VERSION/Organisms/$SPECIES/GeneTranscripts.jdb"))
        catch
            println("Loading GeneTranscripts.jdb failed, replace file with save_RefSeq()")
        end
        global TRANSCRIPTGENE = try
            Dict(load("$PATH/$VERSION/Organisms/$SPECIES/TranscriptGene.jdb"))
        catch
            println("Loading TranscriptGene.jdb failed, replace file with save_RefSeq()")
        end
        global TRANSCRIPTDATA = try
            load("$PATH/$VERSION/Organisms/$SPECIES/TranscriptData.jdb")
        catch
            println("Loading TranscriptData.jdb failed.")
        end
        Base.GC.enable(true);
        println("Loaded $SPECIES RefSeq Successfully")
    else
        println("Species not properly loaded")
    end
end

function Update_Version(version::String)
    write("$(PATH)/RefSeq_version.txt", version);
end

function Check_NCBI_Version()::String
    replace(String(HTTP.get("ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER").body), "\n" => "")
end

function download_RefSeq(version::String=VERSION)
    if version == Check_NCBI_Version()
        println("Current version is most up to date")
    else
        Update_Version(Check_NCBI_Version())
        Load_Version()
        !(isdir("$PATH/$VERSION")) && mkdir("$PATH/$VERSION")
        !(isdir("$PATH/$VERSION/Download")) && mkdir("$PATH/$VERSION/Download")
        link = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.x.rna.gbff.gz"
        i = 1
        iMax = 282
        while true
            try
                println("Downloading vertebrate_mammalian.$(i).rna.gbff.gz")
                download(replace(link, "x" => i), "$PATH/$VERSION/Download/vertebrate_mammalian.$(i).rna.gbff.gz")
            catch
                (i > iMax) && break
            end
            i = i + 1
        end
    end
end

function initial_process_RefSeq()::Int
    counter = 0
    i = 0
    tbl = table((Organism=String[], Name=String[], Transcript=String[], Range=UnitRange{Int64}[], Type=String[], Gene=String[], GeneID=Int[], Description=String[], RefSeq=ReferenceSequence[]))
    !(isdir("$PATH/$VERSION/Processed")) && mkdir("$PATH/$VERSION/Processed")
    for file in readdir("$PATH/$VERSION/Download")
        f = gzopen("$PATH/$VERSION/Download/$file")
        p = ProgressUnknown("Processing $file ... ")
        while(!eof(f))
            x = readuntil(f, "LOCUS   ")
            if x != ""
                organism="";name="";version="";r=1:1;tp="";gene="";id=1; description="";seq=""
                try organism = match(r"ORGANISM  ([^\n]*)", x)[1] catch end
                try name = match(r"SOURCE.*\((.*)\)", x)[1] catch end
                try version = match(r"VERSION[\s]*([NX][RM]\_[\d]*\.[\d])", x)[1] catch end
                try
                    typerange = match(r"((\bCDS\b)|(\bncRNA\b)|(\bmisc_RNA\b)|(\brRNA\b))[\s]*(([\<\>\d]*\.\.[\<\>\d]*)|(join\([\<\>\d]*\.\.[\<\>\d]*\,[\<\>\d]*\.\.[\<\>\d]*))", x)
                    rg = typerange[6]
                    tp = typerange[1]
                    rg = replace(replace(replace(replace(rg, ">" => ""), "<" => ""), "," => ".."), "join(" => "")
                    rgs = split(rg, "..")
                    rg1 = parse(Int, rgs[1])
                    rg2 = parse(Int, rgs[end])
                    r = rg1:rg2
                catch
                end
                try seq = replace(replace(replace(replace(uppercase(match(r"ORIGIN([\s\S]*)\/\/", x)[1]), " " => ""), "\n" => ""), r"[\d]*" => ""), "T" => "U")  catch end
                try gene = match(r"\/gene=\"(\S*)\"", x)[1]  catch end
                try description = match(r"DEFINITION  ([^\(\n]*)", x)[1] catch end
                try id = parse(Int, match(r"\"GeneID:([\d]*)\"", x)[1]) catch end
                push!(rows(tbl), (Organism=organism, Name=name, Transcript=version, Range=r, Type=tp, Gene=gene, GeneID=id, Description=description, RefSeq=encode_refseq(seq)))
                counter=counter+1
                if counter % 100000 == 0
                    i = (Int(trunc(counter / 100000)))
                    Base.GC.enable(false)
                    save(tbl, "$PATH/$VERSION/Processed/DataTbl$i.jdb")
                    Base.GC.enable(true)
                    tbl = table((Organism=String[], Name=String[], Transcript=String[], Range=UnitRange{Int64}[], Type=String[], Gene=String[], GeneID=Int[], Description=String[], RefSeq=ReferenceSequence[]))
                end
            end
            ProgressMeter.next!(p)
        end
        ProgressMeter.finish!(p)
    end
    i = i + 1
    Base.GC.enable(false)
    save(tbl, "$PATH/$VERSION/Processed/DataTbl$i.jdb")
    Base.GC.enable(true)
    return i
end

function secondary_process_RefSeq(i::Int)
    !(isdir("$PATH/$VERSION/Organisms")) && mkdir("$PATH/$VERSION/Organisms")
    p = ProgressMeter.Progress(i, 0.1, "Processing ... ")
    for j in 1:i
        Base.GC.enable(false)
        cur_tbl = load("$PATH/$VERSION/Processed/DataTbl$j.jdb")
        Base.GC.enable(true)
        organisms = JuliaDB.select(cur_tbl, :Organism) |> unique
        r = [replace(x, ".rdb" => "") for x in readdir("$PATH/$VERSION/Organisms")]
        for o in organisms
            println(o)
            o_tbl = filter(r -> r.Organism == o, cur_tbl)
            if o in r
                prev_tbl = MemPool.deserialize("$PATH/$VERSION/Organisms/$(o).rdb")
                o_tbl = merge(o_tbl, prev_tbl)
                prev_tbl = nothing
            end
            MemPool.serialize("$PATH/$VERSION/Organisms/$(o).rdb", o_tbl)
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
end

function tertiary_process_RefSeq()
    for file in readdir("$PATH/$VERSION/Organisms")
        (file[end-3:end] != ".rdb") && continue
        data = MemPool.deserialize("$PATH/$VERSION/Organisms/$(file)")
        !(isdir("$PATH/$VERSION/Organisms/$(replace(file, ".rdb" => ""))")) && mkdir("$PATH/$VERSION/Organisms/$(replace(file, ".rdb" => ""))")
        TranscriptGene = JuliaDB.select(data, (:Transcript, :Gene))
        GeneTranscripts = Dict{String, Array{String, 1}}()
        for (transcript, gene) in TranscriptGene
            if !(haskey(GeneTranscripts, gene))
                GeneTranscripts[gene] = [transcript]
            else
                push!(GeneTranscripts[gene], transcript)
            end
        end
        allRefSeq = JuliaDB.select(data, (:Transcript, :RefSeq))
        TranscriptData = JuliaDB.select(data, (:Gene, :GeneID, :Description, :Transcript, :Range, :Type))
        MemPool.serialize("$PATH/$VERSION/Organisms/$(replace(file, ".rdb" => ""))/TranscriptGene.jdb", TranscriptGene)
        MemPool.serialize("$PATH/$VERSION/Organisms/$(replace(file, ".rdb" => ""))/GeneTranscripts.jdb", GeneTranscripts)
        MemPool.serialize("$PATH/$VERSION/Organisms/$(replace(file, ".rdb" => ""))/allRefSeq.jdb", allRefSeq)
        MemPool.serialize("$PATH/$VERSION/Organisms/$(replace(file, ".rdb" => ""))/TranscriptData.jdb", TranscriptData)
        mv("$(PATH)/$VERSION/Organisms/$(file)", "$(PATH)/$VERSION/Organisms/$(replace(file, ".rdb" => ""))/$(file)")
    end
end
"""
    download_RefSeq(::UnitRange{Int64}=1:8, ::String=PATH)

Downloads mRNA reference sequence from ftp://ftp.ncbi.nlm.nih.gov/refseq/H\\_sapiens/mRNA\\_Prot/ to the PATH folder.  Defaults to downloading 8 files.
"""
function download_RefSeq_old(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
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
function process_RefSeq_old(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
    tbl = table((Transcript=String[], Range=UnitRange{Int64}[], Type=String[], Gene=String[], GeneID=Int[], Description=String[], RefSeq=ReferenceSequence[]))
    for j in num
        f = gzopen("$(path)/$(SPECIES)/$(replace(FILEFORMAT, "X" => "$j"))")
        p = ProgressUnknown("Processing $(replace(FILEFORMAT, "X" => "$j")) ... ")
        while(!eof(f))
            x = readuntil(f, "LOCUS   ")
            try
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
            catch
                println(j)
            end
        end
        ProgressMeter.finish!(p)
        close(f)
    end
    Base.GC.enable(false)
    save(tbl, "$(path)/$(SPECIES)/$(SPECIES)_mRNA_table.jdb")
    Base.GC.enable(true)
    println("Processed data saved to $(path)/$(SPECIES)/$(SPECIES)_mRNA_table")
end


function process_RefSeq_old_old(num::UnitRange{Int64} = 1:NUM, path::String=PATH)
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
function save_RefSeq_old(path::String=PATH)
    println("This may take some time")
    println("--Loading processed data from JDB--")
    Base.GC.enable(false)
    t = load("$(path)/$(SPECIES)/$(SPECIES)_mRNA_table.jdb")
    Base.GC.enable(true)
    println("--Saving Transcript -> Gene Dictionary--")
    TranscriptGene = JuliaDB.select(t, (:Transcript, :Gene))
    Base.GC.enable(false)
    save(TranscriptGene,"$(path)/$(SPECIES)/$(SPECIES)_mRNA_TranscriptGene.jdb")
    Base.GC.enable(true)
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
