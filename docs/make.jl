using Documenter
using siRNATools

makedocs(
    sitename = "siRNATools",
    format = Documenter.HTML(),
    modules = [siRNATools]
)

deploydocs(
    repo = "github.com/maxpwilson/siRNATools.jl.git"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() inthe Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
