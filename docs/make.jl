using Documenter
#using siRNATools

makedocs(
    sitename = "siRNATools",
    format = Documenter.HTML(),
    modules = [siRNATools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
