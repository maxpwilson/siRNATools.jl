
function rnafold(pat::String, file::String="$PATH/rnafold_temp.txt")::Tuple{String, Float64}
    open(file, "w") do io
        write(io, pat)
    end
    raw_output::String = read(`RNAfold $file`, String)
    folding_info::Vector{String} = split(split(raw_output, "\n")[2], " ")
    folding_sites::String = folding_info[1]
    free_energy::Float64 = parse(Float64, replace(folding_info[end], r"\(|\)"=>""))
    (folding_sites, free_energy)
end

function rnafold_firstbond(folding_sites::String)::Int
    firstbond = findfirst("(", folding_sites)
    typeof(firstbond) != Nothing ? firstbond[1] : 0
end

function rnafold_numbonds(folding_sites::String)::Int
    length(folding_sites) - length(replace(folding_sites, "("=>""))
end

function rnafold_bonds(folding_sites::String)::Vector{String}
    bonds::Vector{String} = []
    first_base::Vector{Int} = []
    for x in 1:length(folding_sites)
        if folding_sites[x] == '('
            push!(first_base, x)
        elseif folding_sites[x] == ')'
            push!(bonds, "$(pop!(first_base))-$x")
        end
    end
    bonds
end

function batch_rnafold(pats)::DataFrame
    df = DataFrame(Pattern=String[],Folding=String[],Free_Energy=Float64[],First_Bond=Int[],Num_Bonds=Int[],All_Bonds=String[])
    for pat in pats
        pat_ss = secondary_structure(pat)
        push!(df, [pat, pat_ss.folding_sites,pat_ss.free_energy,pat_ss.first_bond,pat_ss.num_bonds, join(pat_ss.bonds,",")])
    end
    df
end
