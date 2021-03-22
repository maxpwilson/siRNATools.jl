
function Load_SNP_Version()
    try
        if !("SNP_version.txt" in readdir(PATH))
            touch("$(PATH)/SNP_version.txt")
        end
        global SNP_VERSION = readline(open("$(PATH)/SNP_version.txt"))
    catch
    end
end




function process_SNP_db(file::String)
    g = gzopen(file)
    r = readline(g)
    tbl = table((Chrom=String[], Pos=String[], ID=String[], Ref=String[], Alt=String[], Qual=String[], Filter=String[], VC=String[], Freq=String[], FreqYN=String[], Freq1=String[]))
    while r[1:1] == "#"
        r = readline(g)
    end
    counter = 0
    i = 1
    while !eof(g)
        s = split(replace(r, "\n" => ""), "\t")
        chrom = s[1]
        pos = s[2]
        id = s[3]
        ref = s[4]
        alt = s[5]
        qual = s[6]
        filter = s[7]
        vc = !(typeof(s[8]) == Nothing) ? (occursin("VC", s[8]) ? match(r"VC=([^;]*)", s[8])[1] : "") : ""
        freq = !(typeof(s[8]) == Nothing) ? (occursin("FREQ", s[8]) ? match(r"FREQ=([^;]*)", s[8])[1] : "") : ""
        freqyn = (freq == "") ? "No" : "Yes"
        freq1 = "No"
        if freqyn == "Yes"
            for cur in split(freq, "|")
                nums = split(cur, ":")[2]
                n1 = 0
                n2 = 0
                try
                    n1 = parse(Float64, split(nums, ",")[1])
                    n2 = parse(Float64, split(nums, ",")[2])
                catch
                end
                if (n1 >= 0.5 && n2 >= 0.01) || (n2 >= 0.5 && n1 >= 0.01)
                    freq1 = "Yes"
                    break
                end
            end
        end
        push!(rows(tbl), (Chrom=chrom, Pos=pos, ID=id, Ref=ref, Alt=alt, Qual=qual, Filter=filter, VC=vc, Freq=freq, FreqYN=freqyn, Freq1=freq1))
        r = readline(g)
        counter += 1
        (counter % 1000000 == 0) && (println(counter))
        if counter % 10000000 == 0
            Base.GC.enable(false)
            save(tbl, "$PATH/SNPs/154/Processed/DataTbl$i.jdb")
            tbl = table((Chrom=String[], Pos=String[], ID=String[], Ref=String[], Alt=String[], Qual=String[], Filter=String[], VC=String[], Freq=String[]))
            Base.GC.enable(true)
            i += 1
        end
    end
    Base.GC.enable(false)
    save(tbl, "$PATH/SNPs/154/Processed/DataTbl$i.jdb")
    Base.GC.enable(true)
    nothing
end
