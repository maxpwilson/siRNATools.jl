
function Update_ENSEMBL_Version(version)
    write("$(PATH)/Ensembl/ENSEMBL_version.txt", version)
end

function Check_ENSEMBL_Version()::String
    match(r"Ensembl Release (\d*) Databases", String(HTTP.get("ftp://ftp.ensembl.org/pub/current_README").body))[1]
end

function download_ensembl(species)
    ver = Check_ENSEMBL_Version()
    if ENSEMBL_VERSION != ver
        Update_ENSEMBL_Version(ver)
        mkdir("$(PATH)/Ensembl/$(ver)")
        mkdir("$(PATH)/Ensembl/$(ver)/Download")
        mkdir("$(PATH)/Ensembl/$(ver)/Organisms")
    end
    ftp = FTP("ftp://ftp.ensembl.org/pub/current_fasta/$(species)/cdna/")
    files = readdir(ftp)
    for file in files
        (occursin("all.fa.gz", file)) && (download(ftp, file, "$(PATH)/Ensembl/$(ver)/Download/$(species).all.fa.gz"))
    end
end

function process_ensembl(file, vers, species)
    (ENSEMBL_VERSION != vers) && Update_ENSEMBL_Version(vers)
    f = gzopen("$(PATH)/Ensembl/$(vers)/Download/$(file)")
    r = read(f, String)
    rs = split(r, ">")[2:end]
    tbl = table((Accession=String[], Gene=String[], Description=String[], Sequence=String[]))
    p = ProgressMeter.Progress(length(rs), 0.1, "Processing $file ... ")
    for x in rs
        A = ""; G = ""; D = ""; S = ""
        try
            A = match(r"(\S*)", x)[1]
        catch
        end
        try
            G = match(r"gene_symbol:(\S*)", x)[1]
        catch
        end
        try
            D = match(r"description:([^\[]*)", x)[1]
        catch
        end
        try
            S = replace(replace(match(r"\]\n([^ ]*)", x)[1], "\n" => ""), "T" => "U")
        catch
        end
        push!(rows(tbl), (Accession=A, Gene=G, Description=D, Sequence=S))
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
    close(f)
    MemPool.serialize("$(PATH)/Ensembl/$(vers)/Organisms/$(species).jdb", tbl)
end
