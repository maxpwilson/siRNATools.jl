function separate_chromosomes(file::String)
    f = gzopen("$(PATH)Annotations/$(VERSION)/$(file)")
    while !(eof(f))
        r = readuntil(f, "LOCUS   ")
        if r != ""
            if r[5:6] == "NC"
                num = parse(Int, r[8:13])
                open("$(PATH)/Annotations/$(VERSION)/Chrs/chr$(num).gb", "w") do io
                    print(io, r)
                end
                println("Generated Chromosome $(num) file")
            end
        end
    end
    close(f)
end





function skip_header!(f)
    while !(eof(f))
        r = readline(f)
        if match(r"^FEATURES", r) != nothing
            break
        end
    end
end

function make_df(l, xl)
    df = DataFrame(feature=String[], location=String[])
    for x in unique(l)
        if x != "db_xref" && !(x in names(df))
            insertcols!(df,1, x=>String[])
        end
    end
    for x in unique(xl)
        if !(x in names(df))
            insertcols!(df,1, x=>String[])
        end
    end
    df
end

function gen_lists(file::String)
    f = open("$(PATH)$(file)");
    skip_header!(f)
    l = []
    xl = []
    while !(eof(f))
        r = readline(f)
        m = match(r"\s{21}/(\S*)", r)
        if m != nothing
            mx = split(m[1], "=")
            if mx[1] == "db_xref"
               push!(xl, split(mx[2], ":")[1])
            end
            push!(l, split(m[1], "=")[1])
        end
    end
    close(f)
    xl = replace.(xl, "\"" => "")
    (l, xl)
end

function annotate_chromosome(num, l, xl, bSave=true)
    f = open("$(PATH)/Annotations/$(VERSION)/Chrs/chr$(num).gb");
    skip_header!(f)
    df = make_df(l, xl)
    dftemp = make_df(l, xl)
    counter = 0
    while !(eof(f))
        r = readline(f)
        if match(r"^\s{5}\S", r) != nothing
            if size(dftemp)[1] > 0
                df = vcat(df, dftemp)
            end
            dftemp = make_df(l, xl)
            m = match(r"\s*(\S*)\s*(\S*)", r)
            feature = ""
            location = ""
            if m != nothing
                feature = string(m.captures[1])
                (length(m.captures) > 1) && (location = string(m.captures[2]))
            end
            push!(dftemp, ["" for _ in 1:length(names(dftemp))])
            while location[end:end]==","
                r = readline(f)
                m = match(r"\s*(\S*)", r)
                if m != nothing
                    location = location * string(m[1])
                end
            end
            dftemp.feature = [feature]
            dftemp.location = [location]
        elseif match(r"^\s{21}\S", r) != nothing
            m = match(r"\s*/(\S*)", r)
            if m != nothing
                m = replace(m[1], "\"" => "")
                property = split(m, "=")[1]
                value = occursin("=", m) ? split(m, "=")[2] : "TRUE"
                if property == "db_xref"
                    property = split(value, ":")[1]
                    value = split(value, ":")[end]
                end
                if property in names(df)
                    dftemp[end, property] = value
                else

                end
            end
        end
        counter += 1
        if counter % 500000 == 0
            println(counter)
        end
    end
    close(f)
    df = vcat(df, dftemp)
    for col in names(df)
        try
            df[ismissing.(df[col]), col] = ""
        catch
            try
                df[ismissing.(df[col]), col] = -1
            catch
                try
                    df[ismissing.(df[col]), col] = 0
                catch
                    println("$col")
                end
            end
        end
    end
    (bSave) && (CSV.write("$(PATH)/Annotations/$(VERSION)/Chr_Annotations/chromosome$(num)_annotations.csv", df))
    df
end
function annotate_genome()
    df = annotate_chromosome(1, gen_lists("/Annotations/$(VERSION)/Chrs/chr1.gb")...)
    for num in 2:22
        println("Finished annotation chr $(num-1)")
        df = vcat(df, annotate_chromosome(num, gen_lists("/Annotations/$(VERSION)/Chrs/chr1.gb")...))
    end
    println("Finished annotation chr 24")
    CSV.write("$(PATH)/Annotations/$(VERSION)/genome_annotations.csv", df)
end
function create_chrAll()
    df = CSV.read("$(PATH)/Annotations/$(VERSION)/genome_annotations.csv", DataFrame)
    df_t = df[.!ismissing.(df.transcript_id), :]
    df_out = DataFrame(Type=String[], Gene=String[], Transcript=String[], ChrRange=String[], Chromosome=String[], GeneID=Int64[])
    for (type, gene, transcript, location, chromosome, geneid) in zip(df_t.feature, df_t.gene, df_t.transcript_id, df_t.location, df_t.chromosome, df_t.GeneID)
        m = match(r"(complement\(|)(join\(|)([^\)]*)", location)
        md = m[1] == "" ? 1 : -1
        rgs = ""
        (ismissing(chromosome)) && (chromosome = "NA")
        for x in split(m[3], ",")
            xs = split(x, "..")
            rg1 = ""
            rg2 = ""
            try
                rg1 = "$(xs[1]):$(md):$(xs[2])"
                rg2 = "$(xs[2]):$(md):$(xs[1])"
            catch
                rg1 = x
                rg2 = x
            end
            if md == 1
                (rgs != "") && (rgs = rgs * ",")
                rgs = rgs * rg1
            else
                (rgs != "") && (rgs = "," * rgs)
                rgs = rg2 * rgs
            end
        end
        push!(df_out, [type, gene, transcript, rgs, chromosome, geneid])
    end
    CSV.write("$(PATH)chrAll.csv", df_out)
end


function chr_location(srg, pos)
    m = match(r"(complement\(|)(join\(|)([^\)]*)", srg)
    md = m[1] == "" ? 1 : -1
    rgs = []
    for x in split(m[3], ",")
        xs = split(x, "..")
        rg1 = parse(Int, xs[1]):md:parse(Int, xs[2])
        rg2 = parse(Int, xs[2]):md:parse(Int, xs[1])
        rg = md == 1 ? rg1 : rg2
        push!(rgs, rg)
    end
    (md == -1) && (rgs = reverse(rgs))
    counter = 1
    out_pos = -1
    for rg in rgs
        for n in rg
            if n == pos
                out_pos = counter
            else
                counter += 1
            end
        end
        if out_pos > 0
            break
        else
            counter += 1
        end
    end
    out_pos
end
function specialized_annotations()
    df = CSV.read("$(PATH)/Annotations/$(VERSION)/genome_annotations.csv", DataFrame)
    df_transcripts = df[.!ismissing.(df.transcript_id), :]
    CSV.write("$(PATH)/Annotations/$(VERSION)/transcript_annotations.csv", df_transcripts)
    df_gene = df[df.feature .== "gene", :]
    locs = simple_location_parsing.(df_gene.location)
    df_gene.loc1 = [length(i) > 0 ? i[1] : 0 for i in locs]
    df_gene.loc2 = [length(i) > 0 ? i[end] : 0 for i in locs]
    CSV.write("$(PATH)/Annotations/$(VERSION)/gene_annotations.csv", df_gene)
end

function simple_location_parsing(str::String)
    try
        str = replace(str, r">|<"=>"")
        str = replace(replace(str, "complement(" => ""), ")" => "")
        x1 = parse(Int, split(str, "..")[1])
        x2 = parse(Int, split(str, "..")[2])
        rg = (x1:x2)
        rg
    catch
        println("Location parsing failed '$(str)'")
    end
end

function chrom_acc_to_num(str::String)
    try
        m = match(r"NC_([^\.]*)", str)[1]
        parse(Int, m)
    catch
        0
    end
end
