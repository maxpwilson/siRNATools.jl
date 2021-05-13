
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
            tbl = table((Chrom=String[], Pos=String[], ID=String[], Ref=String[], Alt=String[], Qual=String[], Filter=String[], VC=String[], Freq=String[], FreqYN=String[], Freq1=String[]))
            Base.GC.enable(true)
            i += 1
        end
    end
    Base.GC.enable(false)
    save(tbl, "$PATH/SNPs/154/Processed/DataTbl$i.jdb")
    Base.GC.enable(true)
    i
end

function create_gene_SNP_files()
    df = CSV.read("$(PATH)/Annotations/$(VERSION)/gene_annotations.csv", DataFrame)
    for x in 1:size(df)[1]
        df_gene = DataFrame(:Chrom=>String[], :Pos=>String[], :ID=>String[], :Ref=>String[], :Alt=>String[], :Qual=>String[], :Filter=>String[], :VC=>String[], :Freq=>String[], :FreqYN=>String[], :Freq1=>String[])
        CSV.write("$(PATH)/SNPs/154/Genes/$(df[x, :GeneID]).csv", df_gene)
    end
end

function secondary_process_SNP_db(nums::Int)
    df = CSV.read("$(PATH)/Annotations/$(VERSION)/gene_annotations.csv", DataFrame)
    for i in 1:nums
        d = Dict()
        for x in readdir("$(PATH)/SNPs/154/Genes/")
            d[x] = ""
        end
        Base.GC.enable(false)
        tbl = load("$(PATH)/SNPs/154/Processed/DataTbl$(i).jdb")
        Base.GC.enable(true)
        p = ProgressMeter.Progress(length(tbl), 0.1, "Processing tbl $(i)...")
        for x in tbl
            pos = parse(Int, x.Pos)
            genes = df[(df.chromosome .== chrom_acc_to_num(x.Chrom)) .& (df.loc1 .<= pos) .& (df.loc2 .>= pos), :GeneID]
            for gene in genes
                nm = "$(gene).csv"
                if nm in keys(d)
                    if d[nm] == ""
                        d[nm] = CSV.read("$(PATH)/SNPs/154/Genes/$(nm)", DataFrame; types=[String for _ in 1:11])
                    end
                    try
                        push!(d[nm], x)
                    catch
                        println("Failed in gene $(gene) tbl $(x)")
                    end
                end
            end
            ProgressMeter.next!(p)
        end
        ProgressMeter.finish!(p)
        for (x, df_t) in d
            if df_t != ""
                CSV.write("$(PATH)/SNPs/154/Genes/$(x)", df_t)
            end
        end
    end
end
