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
    xl = replace.(xl, "\"" => "")
    (l, xl)
end

function annotate_chromosome(num, l, xl, bSave=true)
    f = open("$(PATH)/Annotations/Chr_Downloads/chr$(num).gb");
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
            dftemp.feature = feature
            dftemp.location = location
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
    df = vcat(df, dftemp)
    for col in names(df)
        try
            df[ismissing.(df[col]), col] = ""
        catch
            try
                df[ismissing.(df[col]), col] = -1
            catch
                df[ismissing.(df[col]), col] = 0
            end
        end
    end
    (bSave) && (CSV.write("$(PATH)/Annotations/Chr_Annotations/chromosome$(num)_annotations.csv", df))
    df
end
function annotate_genome()
    df = annotate_chromosome(1, gen_lists("/Annotations/Chr_Downloads/chr1.gb")...)
    for num in 2:24
        println("Finished annotation chr $(num-1)")
        df = vcat(df, annotate_chromosome(num, gen_lists("/Annotations/Chr_Downloads/chr1.gb")...))
    end
    println("Finished annotation chr 24")
    CSV.write("$(PATH)/Annotations/genome_annotations.csv", df)
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
