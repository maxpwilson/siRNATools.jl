function find_longname(vern)::String
    df = CSV.read("$(PATH)/$(VERSION)/Organisms/Names.csv", DataFrame)
    d = Dict()
    for x in 1:length(df.Longname)
        if !ismissing(df.Vernacular[x])
            d[df.Vernacular[x]] = df.Longname[x]
            for y in split(df.Vernacular[x], ",")
                d[y] = df.Longname[x]
            end
        end
    end
    if vern in keys(d)
        d[vern]
    else
        "Not found"
    end
end
function find_vernacular(longname) ::String
    df = CSV.read("$(PATH)/$(VERSION)/Organisms/Names.csv", DataFrame)
    out = "Not found"
    o = ""
    if longname in df.Longname
        o = df[df.Longname.==longname, 2][1]
    end
    (!ismissing(o)) && (out = o)
    out
end
function search_longnames(longname)
        df = CSV.read("$(PATH)/$(VERSION)/Organisms/Names.csv", DataFrame)
        out = []
        for x in df.Longname
            if !ismissing(x) && occursin(lowercase(longname), lowercase(x))
                push!(out, x)
            end
        end
        out
end

function search_vernacular(vern)
    df = CSV.read("$(PATH)/$(VERSION)/Organisms/Names.csv", DataFrame)
    out = []
    for x in df.Vernacular
        if !ismissing(x) && occursin(lowercase(vern), lowercase(x))
            push!(out, x)
        end
    end
    out
end
