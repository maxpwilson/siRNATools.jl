
function Load_miRBASE_Version()
    try
        if !("miRBASE_version.txt" in readdir(PATH))
            touch("$(PATH)/miRBASE_version.txt")
        end
        global MIRBASE_VERSION = readline(open("$(PATH)/miRBASE_version.txt"))
    catch
    end
end

function Update_miRBASE_Version(version::String)
    write("$(PATH)/miRBASE_version.txt", version);
end

function Check_miRBASE_Version()::String
    r = readline(open(download("ftp://mirbase.org/pub/mirbase/CURRENT/README")))
    try
        match(r"Release ([0-9\.]*)", r)[1]
    catch
        "NA"
    end
end


function download_miRBASE()
    ver = Check_miRBASE_Version()
    if !(ver in readdir("$(PATH)/miRNA/"))
        mkdir("$(PATH)/miRNA/$ver")
        mkdir("$(PATH)/miRNA/$ver/Download")
        mkdir("$(PATH)/miRNA/$ver/Processed")
    end
    download("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz", "$(PATH)/miRNA/$ver/Download/mature_miRBASE.fa.gz")
    Update_miRBASE_Version(ver)
end
function process_miRBASE()
    f = gzopen("$(PATH)/miRNA/$MIRBASE_VERSION/Download/mature_miRBASE.fa.gz")
    r = read(f, String)
    s = split(r, ">")[2:end]
    df = DataFrame(PrevID=String[], ID=String[], Organism=String[], CommonName=String[], Family=String[], Seq=String[])
    for x in s
        previd = try match(r"(\S*)", x)[1]
        catch
            "NA"
        end
        id = try match(r"\S*\s(\S*)", x)[1]
        catch
            "NA"
        end
        organism = try match(r"([A-Z][a-z]*\s[a-z]*)", x)[1]
        catch
            "NA"
        end
        family = try match(r"([^\n\s]*)\n", x)[1]
        catch
            "NA"
        end
        seq = try match(r"\n([^\n\s]*)\n", x)[1]
        catch
            "NA"
        end
        push!(df, [previd, id, organism, split(find_vernacular(organism), ",")[1], family,seq])
    end

    CSV.write("$(PATH)/miRNA/$(MIRBASE_VERSION)/Processed/miRBASE.csv", df)
end

function load_miRBASE()
    global MIRBASE = try
        CSV.read("$(PATH)/miRNA/$(MIRBASE_VERSION)/Processed/miRBASE.csv", DataFrame)
    catch
        println("miRBASE could not be loaded")
    end
end

Load_miRBASE_Version()
load_miRBASE()
