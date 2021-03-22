

function Check_NCBI_Version()::String
    replace(String(HTTP.get("ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER").body), "\n" => "")
end

function Update_NCBI_Version(version::String)
    write("$(PATH)/RefSeq_version.txt", version);
end


function download_RefSeq(version::String=VERSION)
    if version == Check_NCBI_Version()
        println("Current version is most up to date")
    else
        Update_NCBI_Version(Check_NCBI_Version())
        Load_Version()
        !(isdir("$PATH/$VERSION")) && mkdir("$PATH/$VERSION")
        !(isdir("$PATH/$VERSION/Download")) && mkdir("$PATH/$VERSION/Download")
        link = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/vertebrate_mammalian.x.rna.gbff.gz"
        i = 1
        iMax = 400
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

function process_RefSeq(; primary::Bool=true, secondary::Bool=true, tertiary::Bool=true, i::Int=0)
    (primary == true) && (i = initial_process_RefSeq())
    (secondary == true) && (secondary_process_RefSeq(i))
    (tertiary == true) && (tertiary_process_RefSeq())
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
        close(f)
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
