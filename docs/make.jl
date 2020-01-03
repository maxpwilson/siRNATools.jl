using Documenter
using siRNATools

makedocs(
    sitename = "siRNATools",
    format = Documenter.HTML(),
    modules = [siRNATools],
    pages = [
        "Home" => "index.md",
        "Installation Guide" => "man/instal_guide.md",
        "Specificity" => Any[
            "man/spec_guide.md",
            "man/spec_index.md"
        ],
        "Activity" => Any[
            "man/act_guide.md",
            "man/act_index.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/maxpwilson/siRNATools.jl.git"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
